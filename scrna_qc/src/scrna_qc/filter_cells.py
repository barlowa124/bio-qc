"""Threshold filtering with an explicit waterfall.

Every step records cells/genes in and out so the report can show what the
filters removed, not just what survived. Retained + removed must always
reconcile to the input count.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np

from .qc_metrics import add_qc_metrics
from .util import load_config


def apply_filters(adata: ad.AnnData, cfg: dict) -> tuple[ad.AnnData, dict]:
    q = cfg["qc"]
    adata = add_qc_metrics(adata, q["mt_prefix"])

    waterfall = {"cells_in": int(adata.n_obs), "genes_in": int(adata.n_vars),
                 "steps": []}

    def record(stage, before_obs, before_vars):
        waterfall["steps"].append({
            "stage": stage,
            "cells_before": int(before_obs),
            "cells_after": int(adata.n_obs),
            "cells_removed": int(before_obs - adata.n_obs),
            "genes_before": int(before_vars),
            "genes_after": int(adata.n_vars),
            "genes_removed": int(before_vars - adata.n_vars),
        })

    keep = np.asarray(adata.obs["n_genes_by_counts"]) >= q["min_genes_per_cell"]
    adata = adata[keep].copy()
    record("min_genes_per_cell", len(keep), adata.n_vars)

    keep = np.asarray(adata.obs["pct_mt"]) <= q["max_pct_mt"]
    before_v = adata.n_vars
    adata = adata[keep].copy()
    record("max_pct_mt", len(keep), before_v)

    gene_keep = np.asarray((adata.X > 0).sum(axis=0)).ravel() >= q["min_cells_per_gene"]
    before_c = adata.n_obs
    adata = adata[:, gene_keep].copy()
    waterfall["steps"].append({
        "stage": "min_cells_per_gene",
        "cells_before": int(before_c), "cells_after": int(adata.n_obs),
        "cells_removed": int(before_c - adata.n_obs),
        "genes_before": int(len(gene_keep)), "genes_after": int(adata.n_vars),
        "genes_removed": int(len(gene_keep) - adata.n_vars),
    })

    waterfall["cells_out"] = int(adata.n_obs)
    waterfall["genes_out"] = int(adata.n_vars)
    return adata, waterfall


def main() -> None:
    src, dst, waterfall_out = sys.argv[1:4]
    cfg = load_config()
    adata = ad.read_h5ad(src)
    adata, waterfall = apply_filters(adata, cfg)
    Path(dst).parent.mkdir(parents=True, exist_ok=True)
    Path(waterfall_out).parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(dst)
    with open(waterfall_out, "w") as f:
        json.dump(waterfall, f, indent=2)
    print(f"filter: {waterfall['cells_in']} -> {waterfall['cells_out']} cells, "
          f"{waterfall['genes_in']} -> {waterfall['genes_out']} genes")


if __name__ == "__main__":
    main()
