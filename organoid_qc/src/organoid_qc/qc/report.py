"""Report: joint UMAP (reference + organoid) colored by origin and cluster."""

from __future__ import annotations

import json
import sys

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
from organoid_qc.util import load_config


def report(
    organoid: ad.AnnData,
    reference: ad.AnnData,
    fidelity: dict,
    out_path: str,
    cfg: dict,
) -> None:
    ct_col = cfg["dataset"]["obs_columns"]["cell_type"]

    joint = ad.concat(
        {
            "reference": reference,
            "organoid": organoid,
        },
        label="origin",
        index_unique="-",
    )
    sc.pp.neighbors(joint, use_rep="X_pca")
    sc.tl.umap(joint)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    origin = joint.obs["origin"].astype(str)
    colors = {"reference": "#999999", "organoid": "#d62728"}
    axes[0].scatter(
        joint.obsm["X_umap"][:, 0],
        joint.obsm["X_umap"][:, 1],
        c=[colors[o] for o in origin],
        s=4,
        alpha=0.5,
    )
    axes[0].set_title("origin")
    for ct in sorted(reference.obs[ct_col].unique()):
        m = (origin == "reference") & (joint.obs[ct_col] == ct)
        axes[0].scatter(
            joint.obsm["X_umap"][m, 0],
            joint.obsm["X_umap"][m, 1],
            s=6,
            label=ct,
        )
    axes[0].legend(markerscale=3, fontsize=7)

    unmapped = {
        r["cluster"] for r in fidelity["clusters"] if not r["mapped"]
    }
    is_org = origin == "organoid"
    in_unmapped = is_org & joint.obs["leiden"].isin(unmapped)
    axes[1].scatter(
        joint.obsm["X_umap"][~is_org, 0],
        joint.obsm["X_umap"][~is_org, 1],
        c="#dddddd",
        s=4,
        alpha=0.4,
        label="reference",
    )
    axes[1].scatter(
        joint.obsm["X_umap"][is_org & ~in_unmapped, 0],
        joint.obsm["X_umap"][is_org & ~in_unmapped, 1],
        c="#1f77b4",
        s=5,
        label="mapped organoid",
    )
    axes[1].scatter(
        joint.obsm["X_umap"][in_unmapped, 0],
        joint.obsm["X_umap"][in_unmapped, 1],
        c="#d62728",
        s=5,
        label="unmapped organoid",
    )
    axes[1].set_title(
        f"fidelity: coverage {fidelity['aggregate']['celltype_coverage']}, "
        f"unmapped {fidelity['aggregate']['unmapped_fraction']}"
    )
    axes[1].legend(markerscale=3, fontsize=7)
    for ax in axes:
        ax.set_xticks([])
        ax.set_yticks([])
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    print(f"report: {out_path}")


def main() -> None:
    org_path, ref_path, fid_path, out_path = sys.argv[1:5]
    cfg = load_config()
    with open(fid_path) as f:
        fidelity = json.load(f)
    report(
        ad.read_h5ad(org_path),
        ad.read_h5ad(ref_path),
        fidelity,
        out_path,
        cfg,
    )


if __name__ == "__main__":
    main()
