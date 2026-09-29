"""Standard scRNA-seq preprocessing shared by organoid and reference.

Filter -> normalize -> log1p -> HVG (on the joint gene space) -> PCA ->
Leiden clustering for the organoid. Reference and organoid are processed
with identical parameters so their expression profiles are comparable.

Writes *_pp.h5ad; HVG selection is unioned across both datasets and stored
in var['highly_variable'] so scoring uses one consistent gene space.
"""

from __future__ import annotations

import sys

import anndata as ad
import scanpy as sc
from organoid_qc.util import load_config


def _basic_qc(adata: ad.AnnData, cfg: dict) -> ad.AnnData:
    sc.pp.filter_cells(adata, min_genes=cfg["preprocess"]["min_genes_per_cell"])
    sc.pp.filter_genes(adata, min_cells=cfg["preprocess"]["min_cells_per_gene"])
    sc.pp.normalize_total(adata, target_sum=cfg["preprocess"]["target_sum"])
    sc.pp.log1p(adata)
    return adata


def preprocess(organoid: ad.AnnData, reference: ad.AnnData, cfg: dict):
    organoid = _basic_qc(organoid, cfg)
    reference = _basic_qc(reference, cfg)

    # Shared HVG space: intersect on gene names, union HVG flags.
    shared = organoid.var_names.intersection(reference.var_names)
    if not len(shared):
        raise ValueError(
            "organoid and reference share no gene names; nothing to score"
        )
    organoid = organoid[:, shared].copy()
    reference = reference[:, shared].copy()
    for a in (organoid, reference):
        sc.pp.highly_variable_genes(
            a, n_top_genes=min(cfg["preprocess"]["n_hvg"], a.n_vars)
        )
    hvg = (
        organoid.var["highly_variable"] | reference.var["highly_variable"]
    )
    organoid.var["highly_variable"] = hvg
    reference.var["highly_variable"] = hvg

    for a in (organoid, reference):
        a.layers["lognorm"] = a.X.copy()
        sc.pp.scale(a, max_value=10)
        sc.tl.pca(
            a,
            n_comps=min(cfg["preprocess"]["n_pcs"], a.n_vars - 1, a.n_obs - 1),
            use_highly_variable=True,
        )
        sc.pp.neighbors(a)
        sc.tl.leiden(
            a,
            resolution=cfg["preprocess"]["leiden_resolution"],
            key_added="leiden",
            flavor="igraph",
            n_iterations=2,
            directed=False,
        )
    return organoid, reference


def main() -> None:
    org_in, ref_in, org_out, ref_out = sys.argv[1:5]
    cfg = load_config()
    organoid = ad.read_h5ad(org_in)
    reference = ad.read_h5ad(ref_in)
    organoid, reference = preprocess(organoid, reference, cfg)
    organoid.write_h5ad(org_out)
    reference.write_h5ad(ref_out)
    print(
        f"preprocess: organoid {organoid.n_obs} cells / "
        f"{organoid.obs['leiden'].nunique()} clusters; "
        f"reference {reference.n_obs} cells / "
        f"{reference.obs[cfg['dataset']['obs_columns']['cell_type']].nunique()} types"
    )


if __name__ == "__main__":
    main()
