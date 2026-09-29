import anndata as ad
import numpy as np
import pandas as pd

from organoid_qc.qc.score import score


def _adata(X, obs, genes):
    a = ad.AnnData(
        X=np.asarray(X, dtype=np.float32),
        obs=pd.DataFrame(obs),
        var=pd.DataFrame(index=genes),
    )
    a.layers["lognorm"] = a.X.copy()
    a.var["highly_variable"] = True
    return a


def _toy():
    genes = ["m_a", "m_b", "g1", "g2"]
    cfg = {
        "dataset": {"obs_columns": {"cell_type": "cell_type", "cluster": "leiden"}},
        "markers": {"typeA": ["m_a"], "typeB": ["m_b"]},
        "score": {
            "mapping_correlation_min": 0.3,
            "min_celltype_coverage": 1.0,
            "max_unmapped_fraction": 0.1,
            "min_marker_detection": 0.5,
            "marker_expression_ratio": 0.25,
        },
    }
    # Reference: two well-separated types.
    ref = _adata(
        [
            [5, 0, 1, 1], [5, 0, 1, 0], [5, 0, 0, 1],
            [0, 5, 1, 1], [0, 5, 0, 1], [0, 5, 1, 0],
        ],
        {"cell_type": ["typeA"] * 3 + ["typeB"] * 3},
        genes,
    )
    # Organoid: cluster 0 ~ typeA, cluster 1 ~ typeB, cluster 2 aberrant.
    org = _adata(
        [
            [4, 0, 1, 1], [4, 0, 1, 0],      # cluster 0
            [0, 4, 1, 1], [0, 4, 0, 1],      # cluster 1
            [2, 2, 2, 2], [2, 2, 3, 2],      # cluster 2
        ],
        {"leiden": ["0", "0", "1", "1", "2", "2"]},
        genes,
    )
    return org, ref, cfg


def test_clusters_map_to_correct_types():
    result = score(*_toy())
    by_cluster = {r["cluster"]: r for r in result["clusters"]}
    assert by_cluster["0"]["best_match"] == "typeA"
    assert by_cluster["1"]["best_match"] == "typeB"
    assert by_cluster["0"]["mapped"] and by_cluster["1"]["mapped"]


def test_aberrant_cluster_is_unmapped():
    result = score(*_toy())
    by_cluster = {r["cluster"]: r for r in result["clusters"]}
    assert not by_cluster["2"]["mapped"]
    assert result["aggregate"]["unmapped_fraction"] == round(2 / 6, 4)


def test_coverage_and_flags():
    result = score(*_toy())
    assert result["aggregate"]["celltype_coverage"] == 1.0
    assert result["aggregate"]["per_type_marker_detection"] == {
        "typeA": 1.0,
        "typeB": 1.0,
    }
    # Unmapped fraction 2/6 > 0.1 threshold -> flag fails.
    assert result["qc_flags"]["unmapped_fraction"] is False
    assert result["qc_pass"] is False


def test_missing_reference_type_drops_coverage():
    org, ref, cfg = _toy()
    ref = ref[ref.obs["cell_type"] == "typeA"].copy()
    result = score(org, ref, cfg)
    assert result["aggregate"]["celltype_coverage"] == 1.0
    assert result["reference_types"] == ["typeA"]


def test_empty_marker_panel_fails_qc():
    org, ref, cfg = _toy()
    cfg["markers"] = {}  # no panel -> QC question cannot be asked
    result = score(org, ref, cfg)
    assert result["qc_flags"]["marker_detection"] is False
    assert result["qc_pass"] is False
