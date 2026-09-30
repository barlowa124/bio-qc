"""End-to-end: demo data -> filter -> embed -> report, in a tmp dir.

Uses a smaller config so the suite stays fast; exercises the same
functions the Snakefile entry points call.
"""

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from scrna_qc.data import make_demo
from scrna_qc.embed import embed
from scrna_qc.filter_cells import apply_filters
from scrna_qc.report import cluster_table, report


@pytest.fixture()
def cfg():
    return {
        "dataset": {"mode": "demo", "demo": {
            "seed": 3, "n_cells": 400, "n_genes": 500, "n_mt_genes": 12,
            "type_proportions": [0.5, 0.3, 0.2], "marker_block_size": 20,
            "frac_low_quality": 0.10, "frac_mito_high": 0.08}},
        "qc": {"min_genes_per_cell": 100, "min_cells_per_gene": 3,
               "max_pct_mt": 20.0, "mt_prefix": "MT-"},
        "embed": {"target_sum": 10000, "n_hvg": 250, "n_pcs": 20,
                  "n_neighbors": 12, "leiden_resolution": 0.6,
                  "n_markers": 3, "seed": 3},
    }


@pytest.fixture()
def run(tmp_path, cfg):
    raw = make_demo(cfg)
    filtered, waterfall = apply_filters(raw, cfg)
    embedded, markers = embed(filtered, cfg)
    cluster_csv = str(tmp_path / "cluster_qc.csv")
    summary_json = str(tmp_path / "qc_summary.json")
    umap_png = str(tmp_path / "umap.png")
    summary = report(embedded, markers, waterfall, cfg,
                     cluster_csv, summary_json, umap_png)
    return raw, filtered, embedded, markers, waterfall, summary, \
        Path(umap_png)


def test_filters_remove_injected_cells(run):
    raw, filtered, *_ = run
    assert filtered.n_obs < raw.n_obs
    deg = np.asarray(raw.obs["true_degraded"])
    kept = set(filtered.obs_names)
    removed_low = (deg == "low_genes") & \
        np.array([c not in kept for c in raw.obs_names])
    assert removed_low.sum() > (deg == "low_genes").sum() * 0.7


def test_every_cell_clustered(run):
    _, _, embedded, *_ = run
    assert embedded.obs["leiden"].notna().all()
    assert "X_umap" in embedded.obsm
    assert "X_pca" in embedded.obsm


def test_cluster_table_covers_all(run):
    _, _, embedded, markers, _, summary, _ = run
    table = cluster_table(embedded, markers)
    assert len(table) == embedded.obs["leiden"].nunique()
    assert table["n_cells"].sum() == embedded.n_obs
    assert set(table.columns) == {
        "cluster", "n_cells", "frac_of_total", "median_total_counts",
        "median_n_genes", "median_pct_mt", "top_marker"}


def test_markers_have_scores(run):
    _, _, _, markers, *_ = run
    assert len(markers) > 0
    assert markers["scores"].notna().all()


def test_summary_schema_and_provenance(run):
    *_, waterfall, summary, umap_png = run
    assert summary["n_clusters"] >= 2
    assert summary["filter_waterfall"]["cells_out"] == summary["n_cells"]
    prov = summary["provenance"]
    assert prov["scanpy"] and prov["config"]["dataset"]["mode"] == "demo"
    assert umap_png.exists() and umap_png.stat().st_size > 10000


def test_leiden_deterministic(cfg):
    a = embed(apply_filters(make_demo(cfg), cfg)[0], cfg)[0]
    b = embed(apply_filters(make_demo(cfg), cfg)[0], cfg)[0]
    assert a.obs["leiden"].tolist() == b.obs["leiden"].tolist()
