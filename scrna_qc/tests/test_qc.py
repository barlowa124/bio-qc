"""QC metrics and filter waterfall against a hand-worked fixture."""

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from scrna_qc.filter_cells import apply_filters
from scrna_qc.qc_metrics import add_qc_metrics


def fixture():
    # 4 cells x 5 genes; genes 0-1 are "MT-".
    X = np.array([
        [10, 0, 20, 0, 10],   # 40 counts, 3 genes, 25% mt
        [0, 0, 5, 0, 0],      # 5 counts, 1 gene, 0% mt
        [8, 8, 0, 0, 0],      # 16 counts, 2 genes, 100% mt
        [0, 0, 4, 6, 4],      # 14 counts, 3 genes, 0% mt
    ], dtype=float)
    return ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=[f"c{i}" for i in range(4)]),
        var=pd.DataFrame(index=["MT-1", "MT-2", "G1", "G2", "G3"]),
    )


def test_metrics_exact():
    a = add_qc_metrics(fixture(), "MT-")
    assert a.obs["total_counts"].tolist() == [40, 5, 16, 14]
    assert a.obs["n_genes_by_counts"].tolist() == [3, 1, 2, 3]
    assert np.allclose(a.obs["pct_mt"], [25, 0, 100, 0])


def test_no_mt_genes_gives_zero_pct():
    a = fixture()
    a.var_names = [f"G{i}" for i in range(5)]
    a = add_qc_metrics(a, "MT-")
    assert np.all(a.obs["pct_mt"] == 0)


def cfg_for(**over):
    q = {"min_genes_per_cell": 3, "min_cells_per_gene": 2,
         "max_pct_mt": 50.0, "mt_prefix": "MT-"}
    q.update(over)
    return {"qc": q}


def test_filter_boundaries():
    # min_genes=3 keeps c0 and c4 (3 genes), drops c1 (1) and c2 (2).
    # pct_mt<=50 then keeps both survivors; min_cells=2 drops G2's
    # count check per remaining cells.
    a, wf = apply_filters(fixture(), cfg_for())
    assert a.obs_names.tolist() == ["c0", "c3"]
    assert wf["cells_in"] == 4 and wf["cells_out"] == 2
    step = wf["steps"][0]
    assert step["stage"] == "min_genes_per_cell"
    assert step["cells_removed"] == 2


def test_pct_mt_boundary_inclusive():
    a, _ = apply_filters(fixture(), cfg_for(min_genes_per_cell=1,
                                            max_pct_mt=25.0,
                                            min_cells_per_gene=0))
    # c0 at exactly 25% stays; c2 at 100% is removed.
    assert "c0" in a.obs_names and "c2" not in a.obs_names


def test_waterfall_reconciles():
    a, wf = apply_filters(fixture(), cfg_for())
    cells_removed = sum(s["cells_removed"] for s in wf["steps"])
    genes_removed = sum(s["genes_removed"] for s in wf["steps"])
    assert wf["cells_in"] - cells_removed == wf["cells_out"]
    assert wf["genes_in"] - genes_removed == wf["genes_out"]


def test_empty_result_keeps_schema():
    a, wf = apply_filters(fixture(), cfg_for(min_genes_per_cell=99))
    assert a.n_obs == 0
    assert wf["cells_out"] == 0
    assert "total_counts" in a.obs.columns
