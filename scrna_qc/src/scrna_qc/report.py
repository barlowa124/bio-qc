"""Report: per-cluster QC table, summary JSON with provenance, UMAP PNG.

The cluster table is the unit that matters — a pooled QC pass can hide a
depleted or stressed population, so every cluster gets its own row.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import pandas as pd
import scanpy as sc

from .util import git_sha, load_config


def cluster_table(adata: ad.AnnData, markers: pd.DataFrame) -> pd.DataFrame:
    obs = adata.obs.copy()
    obs["cluster"] = obs["leiden"].astype(str)
    group_col = "group" if "group" in markers.columns else "cluster"
    markers = markers.assign(**{group_col: markers[group_col].astype(str)})
    top_marker = (markers.sort_values("scores", ascending=False)
                  .groupby(group_col, observed=True).head(1)
                  .set_index(group_col)["names"].to_dict())

    rows = []
    for cl, sub in obs.groupby("cluster", observed=True):
        rows.append({
            "cluster": cl,
            "n_cells": int(len(sub)),
            "frac_of_total": round(len(sub) / len(obs), 4),
            "median_total_counts": float(sub["total_counts"].median()),
            "median_n_genes": float(sub["n_genes_by_counts"].median()),
            "median_pct_mt": round(float(sub["pct_mt"].median()), 2),
            "top_marker": top_marker.get(cl),
        })
    return pd.DataFrame(rows).sort_values("cluster").reset_index(drop=True)


def report(adata: ad.AnnData, markers: pd.DataFrame, waterfall: dict,
           cfg: dict, cluster_csv: str, summary_json: str,
           umap_png: str) -> dict:
    table = cluster_table(adata, markers)
    Path(cluster_csv).parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(cluster_csv, index=False)

    summary = {
        "dataset_mode": cfg["dataset"]["mode"],
        "n_cells": int(adata.n_obs),
        "n_clusters": int(adata.obs["leiden"].nunique()),
        "filter_waterfall": waterfall,
        "cluster_qc": table.to_dict(orient="records"),
        "provenance": {
            "git_sha": git_sha(),
            "scanpy": sc.__version__,
            "anndata": ad.__version__,
            "config": cfg,
        },
    }
    with open(summary_json, "w") as f:
        json.dump(summary, f, indent=2)

    sc.pl.umap(adata, color=["leiden"], show=False)
    import matplotlib.pyplot as plt
    plt.gcf().set_size_inches(6, 5)
    plt.gcf().savefig(umap_png, dpi=130, bbox_inches="tight")
    plt.close()
    return summary


def main() -> None:
    emb_in, markers_in, waterfall_in, cluster_csv, summary_json, umap_png = \
        sys.argv[1:7]
    cfg = load_config()
    adata = ad.read_h5ad(emb_in)
    markers = pd.read_csv(markers_in)
    with open(waterfall_in) as f:
        waterfall = json.load(f)
    s = report(adata, markers, waterfall, cfg,
               cluster_csv, summary_json, umap_png)
    print(f"report: {s['n_clusters']} clusters over {s['n_cells']} cells")


if __name__ == "__main__":
    main()
