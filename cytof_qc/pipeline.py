"""End-to-end run: load gated FCS files, QC, transform, cluster, score.

Writes compact artifacts to results/:

- ``metrics.json`` — channel QC, population QC, agreement metrics,
  per-population recovery table
- ``umap_populations.png`` / ``umap_clusters.png``
- ``marker_heatmap.png``

Raw FCS inputs stay out of git; run ``scripts/fetch_data.sh`` first.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import seaborn as sns

from . import cluster, io, mapping, qc, transform

HERE = Path(__file__).resolve().parent.parent
DATA = HERE / "data" / "Levine_13dim"
RESULTS = HERE / "results"


def main(data_dir: Path = DATA) -> dict:
    events = io.load_population_dir(data_dir)
    channels = [c for c in events.columns if c not in {"population", "source_file"}]

    channel_report = qc.channel_qc(events, channels)
    population_report = qc.population_qc(events)

    # cluster every event, gated or not: the benchmark protocol runs the
    # clustering on the full file and scores only the labeled subset
    transformed = transform.arcsinh(events[channels])
    adata = cluster.cluster_events(
        transformed, channels, n_neighbors=15, resolution=1.0
    )
    events["cluster"] = adata.obs["cluster"].to_numpy()
    gated = events[events["population"] != "NotGated"].copy()

    agree = mapping.agreement(
        gated["population"].to_numpy(), gated["cluster"].to_numpy()
    )
    table = mapping.contingency(gated["population"], gated["cluster"])
    cluster_map = mapping.match_clusters_to_populations(table)
    per_pop = mapping.per_population_report(gated["population"], gated["cluster"])

    metrics = {
        "dataset": "Levine_13dim (Levine et al. 2015, human bone marrow, "
        "manually gated populations)",
        "n_events_clustered": int(len(events)),
        "n_events_scored": int(len(gated)),
        "n_channels": len(channels),
        "n_gated_populations": int(gated["population"].nunique()),
        "n_clusters_found": int(events["cluster"].nunique()),
        "agreement": agree,
        "channels": channel_report,
        "populations": population_report,
        "per_population_recovery": per_pop,
        "cluster_to_population": cluster_map,
    }

    RESULTS.mkdir(exist_ok=True)
    (RESULTS / "metrics.json").write_text(json.dumps(metrics, indent=2))

    _figures(adata, events, table)

    summary = (
        f"{metrics['n_events_clustered']} events clustered, "
        f"{metrics['n_events_scored']} scored, "
        f"{metrics['n_gated_populations']} gated populations, "
        f"{metrics['n_clusters_found']} Leiden clusters; "
        f"ARI={agree['adjusted_rand_index']}, NMI={agree['normalized_mutual_info']}"
    )
    print(summary)
    return metrics


def _figures(adata, events: pd.DataFrame, table: pd.DataFrame) -> None:
    adata.obs["population"] = events["population"].to_numpy()

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    sc.pl.umap(adata, color="population", ax=axes[0], show=False, legend_loc=None,
               title="manual gates")
    sc.pl.umap(adata, color="cluster", ax=axes[1], show=False, legend_loc=None,
               title="Leiden clusters")
    fig.savefig(RESULTS / "umap_comparison.png", dpi=120, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 7))
    sns.heatmap(
        table.div(table.sum(axis=1), axis=0), cmap="mako", ax=ax,
        cbar_kws={"label": "fraction of population"},
    )
    ax.set_xlabel("Leiden cluster")
    ax.set_ylabel("gated population")
    fig.savefig(RESULTS / "population_cluster_heatmap.png", dpi=120,
                bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main(Path(sys.argv[1]) if len(sys.argv) > 1 else DATA)
