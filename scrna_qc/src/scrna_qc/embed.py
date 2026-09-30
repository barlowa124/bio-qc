"""Embedding stage: normalize -> log1p -> HVG -> scale -> PCA -> UMAP ->
Leiden -> per-cluster marker detection (Wilcoxon).

All stochastic steps take the config seed so a rerun reproduces the same
cluster assignments.
"""

from __future__ import annotations

import sys
from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc

from .util import load_config


def embed(adata: ad.AnnData, cfg: dict) -> tuple[ad.AnnData, pd.DataFrame]:
    e = cfg["embed"]
    seed = e["seed"]

    sc.pp.normalize_total(adata, target_sum=e["target_sum"])
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=min(e["n_hvg"], adata.n_vars))
    adata.raw = adata
    adata = adata[:, adata.var.highly_variable].copy()
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=min(e["n_pcs"], adata.n_vars - 1), random_state=seed)
    sc.pp.neighbors(adata, n_neighbors=e["n_neighbors"], random_state=seed)
    sc.tl.umap(adata, random_state=seed)
    sc.tl.leiden(adata, resolution=e["leiden_resolution"], random_state=seed)
    sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")

    markers = sc.get.rank_genes_groups_df(adata, group=None)
    group_col = "cluster" if "cluster" in markers.columns else "group"
    markers = markers.groupby(group_col, observed=True).head(e["n_markers"])
    return adata, markers


def main() -> None:
    src, dst, markers_out = sys.argv[1:4]
    cfg = load_config()
    adata = ad.read_h5ad(src)
    adata, markers = embed(adata, cfg)
    Path(dst).parent.mkdir(parents=True, exist_ok=True)
    Path(markers_out).parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(dst)
    markers.to_csv(markers_out, index=False)
    print(f"embed: {adata.n_obs} cells, {adata.obs['leiden'].nunique()} clusters")


if __name__ == "__main__":
    main()
