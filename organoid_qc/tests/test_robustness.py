"""Robustness battery: scoring degenerate profiles, empty gene spaces,
missing markers, and QC flag boundaries."""

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from organoid_qc.automation.simulate import runlog_text
from organoid_qc.qc.score import _corr, _mean_profile, score


def _adata(X, obs, genes, hvg=True):
    a = ad.AnnData(
        X=np.asarray(X, dtype=np.float32),
        obs=pd.DataFrame(obs),
        var=pd.DataFrame(index=genes),
    )
    a.layers["lognorm"] = a.X.copy()
    a.var["highly_variable"] = hvg
    return a


def _cfg(**over):
    base = {
        "dataset": {"obs_columns": {"cell_type": "cell_type",
                                    "cluster": "leiden"}},
        "markers": {"typeA": ["m_a"], "typeB": ["m_b"]},
        "score": {"mapping_correlation_min": 0.3,
                  "min_celltype_coverage": 1.0,
                  "max_unmapped_fraction": 0.1,
                  "min_marker_detection": 0.5,
                  "marker_expression_ratio": 0.25},
    }
    base["score"].update(over)
    return base


def _ref(genes=("m_a", "m_b", "g1", "g2")):
    return _adata(
        [[5, 0, 1, 1], [5, 0, 1, 0], [5, 0, 0, 1],
         [0, 5, 1, 1], [0, 5, 0, 1], [0, 5, 1, 0]],
        {"cell_type": ["typeA"] * 3 + ["typeB"] * 3},
        list(genes))


class TestCorrEdges:
    def test_constant_profile_corr_zero(self):
        a = np.ones(5)
        b = np.arange(5, dtype=float)
        assert _corr(a, b) == 0.0
        assert _corr(a, a) == 0.0  # no 1.0 false confidence

    def test_empty_profiles(self):
        a, b = np.array([]), np.array([])
        with np.errstate(all="ignore"):
            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                v = _corr(a, b)
        assert v == 0.0 or np.isnan(v)  # pin: no crash


class TestMeanProfileEdges:
    def test_missing_genes_raise(self):
        a = _adata([[1, 2]], {"cell_type": ["t"]}, ["g1", "g2"])
        with pytest.raises(ValueError, match="absent"):
            _mean_profile(a, np.array([True]), ["g1", "missing_gene"])

    def test_empty_gene_list_returns_empty(self):
        a = _adata([[1, 2]], {"cell_type": ["t"]}, ["g1", "g2"])
        out = _mean_profile(a, np.array([True]), [])
        assert out.shape == (0,)


class TestScoreEdges:
    def test_no_hvg_flagged_gives_unmapped_not_crash(self):
        org = _adata([[1, 0, 0, 0]] * 4, {"leiden": ["0"] * 4},
                     ["m_a", "m_b", "g1", "g2"], hvg=False)
        ref = _ref()
        ref.var["highly_variable"] = False
        with np.errstate(all="ignore"):
            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                res = score(org, ref, _cfg())
        # empty profile correlations are degenerate -> nothing maps
        assert res["aggregate"]["celltype_coverage"] == 0.0
        assert res["qc_pass"] is False

    def test_identical_constant_cells(self):
        # every cell identical -> all correlations 0 -> all unmapped
        org = _adata(np.ones((6, 4)), {"leiden": ["0", "0", "0",
                                                  "1", "1", "1"]},
                     ["m_a", "m_b", "g1", "g2"])
        res = score(org, _ref(), _cfg())
        assert all(not r["mapped"] for r in res["clusters"])
        assert res["qc_pass"] is False

    def test_unmapped_fraction_boundary(self):
        # one mapped + one unmapped cluster of equal size -> 0.5 unmapped
        org = _adata(
            [[4, 0, 1, 1], [4, 0, 1, 0], [0, 0, 9, 9], [0, 0, 9, 8]],
            {"leiden": ["0", "0", "1", "1"]},
            ["m_a", "m_b", "g1", "g2"])
        res = score(org, _ref(), _cfg(max_unmapped_fraction=0.5))
        assert res["aggregate"]["unmapped_fraction"] == pytest.approx(0.5)
        assert res["qc_flags"]["unmapped_fraction"] is True

    def test_missing_marker_genes_yield_none_not_false(self):
        org = _adata([[4, 0, 1, 1]] * 3, {"leiden": ["0"] * 3},
                     ["m_a", "m_b", "g1", "g2"])
        cfg = _cfg()
        cfg["markers"] = {"typeA": ["no_such_gene"],
                          "typeB": ["also_missing"]}
        res = score(org, _ref(), cfg)
        md = res["aggregate"]["per_type_marker_detection"]
        assert all(v is None for v in md.values())
        # no assessable markers -> flag False (not vacuous pass)
        assert res["qc_flags"]["marker_detection"] is False

    def test_single_cluster_organoid(self):
        org = _adata([[4, 0, 1, 1], [4, 0, 0, 1]],
                     {"leiden": ["0", "0"]},
                     ["m_a", "m_b", "g1", "g2"])
        res = score(org, _ref(), _cfg())
        assert res["aggregate"]["celltype_coverage"] == pytest.approx(0.5)


class TestAutomationEdges:
    def test_empty_runlog(self):
        assert runlog_text([]) == ""

    def test_missing_payload_keys(self):
        txt = runlog_text([{"level": "warning"},
                           {"payload": {"text": "moved 10ul"}}])
        assert "warning" in txt and "moved 10ul" in txt
