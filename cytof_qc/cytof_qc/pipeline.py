"""End-to-end run: load gated FCS files, QC, transform, cluster, score.

Writes compact artifacts per panel under results/<panel>/:

- ``metrics.json`` — channel QC, population QC, drift/batch report,
  agreement metrics, per-population recovery table
- ``umap_populations.png`` / ``umap_clusters.png``
- ``marker_heatmap.png``

Raw FCS inputs stay out of git; run ``scripts/fetch_data.sh`` first.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import seaborn as sns

from . import batch, cluster, io, mapping, qc, transform

HERE = Path(__file__).resolve().parent.parent
DATA = HERE / "data"
RESULTS = HERE / "results"

PANELS = {
    "levine_13dim": {
        "data": "Levine_13dim",
        "label": "Levine_13dim (Levine et al. 2015, human bone marrow, "
        "manually gated populations)",
    },
    "levine_32dim": {
        "data": "Levine_32dim",
        "label": "Levine_32dim (Levine et al. 2015, human bone marrow, "
        "32 markers, two donors H1/H2 acquired on different days)",
    },
}

_META_COLS = {"population", "source_file", "batch"}


def main(
    data_dir: Path,
    out_dir: Path,
    dataset_label: str,
) -> dict:
    events = io.load_population_dir(data_dir)
    all_channels = [c for c in events.columns if c not in _META_COLS]
    channels, qc_only = qc.split_channels(all_channels)

    channel_report = qc.channel_qc(events, all_channels)
    population_report = qc.population_qc(events)

    transformed = transform.arcsinh(events[all_channels])
    transformed["batch"] = events["batch"]

    # acquisition-order drift gate, measured per batch so a batch seam
    # is not mistaken for in-run drift
    has_batch = transformed["batch"].notna().any()
    drift_report = qc.acquisition_qc(
        transformed, all_channels, batch_col="batch" if has_batch else None
    )
    # On gated-population benchmark files the row order follows file order
    # (per-population FCS), so drift here reflects population composition
    # changes between files, not instrument decay. The gate is meaningful
    # for continuous acquisitions; the flags below are reported with that
    # caveat rather than treated as instrument failures.
    drift_scope = (
        "row order follows per-population file order in this benchmark; "
        "drift flags reflect file-to-file composition shifts, not "
        "instrument stability"
    )

    # when files carry a batch suffix (Levine_32dim H1/H2), measure and
    # remove the per-channel between-batch median shift before clustering
    batch_report = None
    if has_batch:
        transformed, batch_report = batch.align_batches(
            transformed, channels
        )

    # cluster every event on the phenotypic markers only; DNA/viability/
    # acquisition channels are QC inputs, not phenotype
    adata = cluster.cluster_events(
        transformed[channels], channels, n_neighbors=15, resolution=1.0
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
        "dataset": dataset_label,
        "n_events_clustered": int(len(events)),
        "n_events_scored": int(len(gated)),
        "n_channels": len(channels),
        "n_gated_populations": int(gated["population"].nunique()),
        "n_clusters_found": int(events["cluster"].nunique()),
        "agreement": agree,
        "channels": channel_report,
        "phenotypic_channels": channels,
        "qc_only_channels": qc_only,
        "acquisition_drift": drift_report,
        "acquisition_drift_scope": drift_scope,
        "populations": population_report,
        "per_population_recovery": per_pop,
        "cluster_to_population": cluster_map,
    }

    # per-batch agreement where replicates exist, so batch mixing is
    # measured rather than assumed away
    if batch_report is not None:
        metrics["batch_alignment"] = {
            ch: shifts for ch, shifts in batch_report.items()
        }
        per_batch = {}
        for b, sub in gated.groupby("batch"):
            per_batch[str(b)] = mapping.agreement(
                sub["population"].to_numpy(), sub["cluster"].to_numpy()
            )
        metrics["agreement_by_batch"] = per_batch

    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "metrics.json").write_text(json.dumps(metrics, indent=2))

    _figures(adata, events, table, out_dir)

    summary = (
        f"{metrics['n_events_clustered']} events clustered, "
        f"{metrics['n_events_scored']} scored, "
        f"{metrics['n_gated_populations']} gated populations, "
        f"{metrics['n_clusters_found']} Leiden clusters; "
        f"ARI={agree['adjusted_rand_index']}, NMI={agree['normalized_mutual_info']}"
    )
    print(summary)
    return metrics


def _figures(adata, events: pd.DataFrame, table: pd.DataFrame, out_dir: Path) -> None:
    adata.obs["population"] = events["population"].to_numpy()

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    sc.pl.umap(adata, color="population", ax=axes[0], show=False, legend_loc=None,
               title="manual gates")
    sc.pl.umap(adata, color="cluster", ax=axes[1], show=False, legend_loc=None,
               title="Leiden clusters")
    fig.savefig(out_dir / "umap_comparison.png", dpi=120, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 7))
    sns.heatmap(
        table.div(table.sum(axis=1), axis=0), cmap="mako", ax=ax,
        cbar_kws={"label": "fraction of population"},
    )
    ax.set_xlabel("Leiden cluster")
    ax.set_ylabel("gated population")
    fig.savefig(out_dir / "population_cluster_heatmap.png", dpi=120,
                bbox_inches="tight")
    plt.close(fig)


def _cli() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--panel", choices=list(PANELS), default=None,
        help="which gated benchmark panel to run",
    )
    ap.add_argument("--data-dir", type=Path, default=None,
                    help="directory of gated FCS files (overrides --panel)")
    ap.add_argument("--out", type=Path, default=None,
                    help="output directory (required with --data-dir)")
    ap.add_argument("--label", default=None, help="dataset label")
    args = ap.parse_args()
    if args.data_dir:
        if not args.out:
            ap.error("--out is required with --data-dir")
        main(data_dir=args.data_dir, out_dir=args.out,
             dataset_label=args.label or args.data_dir.name)
        return
    name = args.panel or "levine_13dim"
    spec = PANELS[name]
    main(
        data_dir=DATA / spec["data"],
        out_dir=RESULTS / name,
        dataset_label=spec["label"],
    )


if __name__ == "__main__":
    _cli()
