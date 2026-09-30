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

from .claims import flatten_results, verify_markdown
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


def _conformal_lines(conf: dict) -> list[str]:
    lines = [
        "",
        "## Cluster-assignment confidence (Mondrian conformal)",
        "",
        f"- alpha: {conf['alpha']}, calibration cells: {conf['n_cal']}, "
        f"held-out cells: {conf['n_test']}",
        f"- pooled prediction-set coverage: {conf['pooled_coverage']}",
        f"- mean set size: {conf['mean_set_size']}",
        "",
        "| cluster | held-out cells | coverage | qhat | set size | "
        "fallback |",
        "|---|---|---|---|---|---|",
    ]
    for cl, row in sorted(conf["per_cluster"].items()):
        lines.append(
            f"| {cl} | {row['n_test']} | {row['coverage']} "
            f"| {row['qhat']} | {row['mean_set_size']} "
            f"| {row['global_fallback']} |")
    return lines


def render_markdown(summary: dict, conformal: dict | None = None) -> str:
    """The human-readable report. Every number here must bind to a
    recorded value — verify_markdown enforces that before this text
    is presented as output."""
    wf = summary["filter_waterfall"]
    lines = [
        f"# scrna_qc report — {summary['dataset_mode']}",
        "",
        f"- cells in: {wf['cells_in']}, genes in: {wf['genes_in']}",
        f"- cells after filtering: {summary['n_cells']}, "
        f"genes after filtering: {wf['genes_out']}",
        "- filter removals (cells / genes): "
        + ", ".join(f"{s['stage']} {s['cells_removed']}/{s['genes_removed']}"
                    for s in wf["steps"]),
        f"- Leiden clusters: {summary['n_clusters']}",
        "",
        "## Per-cluster QC",
        "",
        "| cluster | cells | % of cells | median counts | median genes |"
        " median %mt | top marker |",
        "|---|---|---|---|---|---|---|",
    ]
    for row in summary["cluster_qc"]:
        lines.append(
            f"| {row['cluster']} | {row['n_cells']} "
            f"| {round(row['frac_of_total'] * 100, 2)}% "
            f"| {row['median_total_counts']} | {row['median_n_genes']} "
            f"| {row['median_pct_mt']} | {row['top_marker']} |")
    if conformal:
        lines += _conformal_lines(conformal)
    prov = summary["provenance"]
    lines += [
        "",
        "## Provenance",
        "",
        f"git {prov['git_sha']} · scanpy {prov['scanpy']} · "
        f"anndata {prov['anndata']}",
        "",
    ]
    return "\n".join(lines)


def report(adata: ad.AnnData, markers: pd.DataFrame, waterfall: dict,
           cfg: dict, cluster_csv: str, summary_json: str,
           umap_png: str, report_md: str | None = None,
           claims_json: str | None = None,
           conformal_json: str | None = None) -> dict:
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

    if report_md and claims_json:
        # Claim check: re-derive every number in the markdown from the
        # recorded waterfall+summary(+conformal) values; flag what
        # doesn't bind.
        conformal = None
        if conformal_json and Path(conformal_json).exists():
            conformal = json.loads(Path(conformal_json).read_text())
        md = render_markdown(summary, conformal)
        Path(report_md).write_text(md)
        flat = {**flatten_results(waterfall), **flatten_results(summary)}
        if conformal:
            flat.update(flatten_results(conformal))
        labels = {str(r["cluster"]) for r in table.to_dict("records")}
        labels.add(str(summary["provenance"]["git_sha"]))
        verdict = verify_markdown(md, flat, labels)
        with open(claims_json, "w") as f:
            json.dump(verdict, f, indent=2)
        status = "verified" if verdict["passed"] else "UNBOUND CLAIMS"
        print(f"claims: {verdict['n_claims']} numbers checked — {status}")

    return summary


def main() -> None:
    emb_in, markers_in, waterfall_in, cluster_csv, summary_json, umap_png = \
        sys.argv[1:7]
    report_md = sys.argv[7] if len(sys.argv) > 7 else None
    claims_json = sys.argv[8] if len(sys.argv) > 8 else None
    conformal_json = sys.argv[9] if len(sys.argv) > 9 else None
    cfg = load_config()
    adata = ad.read_h5ad(emb_in)
    markers = pd.read_csv(markers_in)
    with open(waterfall_in) as f:
        waterfall = json.load(f)
    s = report(adata, markers, waterfall, cfg,
               cluster_csv, summary_json, umap_png, report_md, claims_json,
               conformal_json)
    print(f"report: {s['n_clusters']} clusters over {s['n_cells']} cells")


if __name__ == "__main__":
    main()
