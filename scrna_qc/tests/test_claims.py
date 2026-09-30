"""Verified-claims layer: numbers in the markdown report must trace to
recorded JSON values."""

import json

import pandas as pd
import pytest

from scrna_qc.claims import (extract_numbers, flatten_results,
                             verify_markdown)
from scrna_qc.report import render_markdown, report


def _summary():
    return {
        "dataset_mode": "demo",
        "n_cells": 783,
        "n_clusters": 4,
        "filter_waterfall": {
            "cells_in": 900, "genes_in": 1200,
            "steps": [
                {"stage": "min_genes_per_cell", "cells_before": 900,
                 "cells_after": 828, "cells_removed": 72,
                 "genes_before": 1200, "genes_after": 1200,
                 "genes_removed": 0},
                {"stage": "max_pct_mt", "cells_before": 828,
                 "cells_after": 783, "cells_removed": 45,
                 "genes_before": 1200, "genes_after": 1200,
                 "genes_removed": 0},
            ],
            "cells_out": 783, "genes_out": 1200,
        },
        "cluster_qc": [
            {"cluster": "0", "n_cells": 352, "frac_of_total": 0.4496,
             "median_total_counts": 2450.0, "median_n_genes": 812.0,
             "median_pct_mt": 4.82, "top_marker": "g_100"},
            {"cluster": "1", "n_cells": 196, "frac_of_total": 0.2503,
             "median_total_counts": 1980.5, "median_n_genes": 655.0,
             "median_pct_mt": 3.91, "top_marker": "g_350"},
        ],
        "provenance": {"git_sha": "abc1234", "scanpy": "1.11.5",
                       "anndata": "0.11.4", "config": {}},
    }


def _flat(summary):
    return {**flatten_results(summary["filter_waterfall"]),
            **flatten_results(summary)}


# --- flatten_results ---

def test_flatten_skips_bools_and_strings():
    flat = flatten_results({"a": 1, "b": True, "c": "x",
                            "d": {"e": [2.5, 3]}})
    assert flat == {"a": 1.0, "d.e.0": 2.5, "d.e.1": 3.0}


# --- extract_numbers ---

def test_extract_numbers_boundaries():
    toks = extract_numbers("cells 900, 48.65% retained, id x2.5, v1.11.5")
    raws = [t["token"] for t in toks]
    assert "900" in raws and "48.65%" in raws
    # version strings and identifier suffixes are not measurements
    assert not any("1.11" in r or "2.5" == r for r in raws)


def test_extract_numbers_signed_range():
    toks = extract_numbers("range 0.5-0.6 and -0.544 logfc")
    vals = [t["value"] for t in toks]
    assert 0.5 in vals and 0.6 in vals and -0.544 in vals


# --- verify_markdown ---

def test_generated_report_is_fully_bound():
    summary = _summary()
    md = render_markdown(summary)
    labels = {r["cluster"] for r in summary["cluster_qc"]}
    labels.add(summary["provenance"]["git_sha"])
    v = verify_markdown(md, _flat(summary), labels)
    assert v["passed"], v
    assert v["n_claims"] > 20


def test_percent_token_binds_to_fraction():
    flat = {"cluster_qc.0.frac_of_total": 0.4496}
    v = verify_markdown("cluster 0 holds 44.96% of cells", flat,
                        labels={"0"})
    assert v["passed"], v


def test_unbound_number_flagged():
    flat = _flat(_summary())
    v = verify_markdown("we kept 9999 cells", flat)
    assert not v["passed"]
    assert v["unbound_claims"][0]["token"] == "9999"


def test_sci_notation_tokens_are_single_claims():
    # "1.2e-20" must not fragment into "1" and "20"
    toks = extract_numbers("p = 1.2e-20")
    assert len(toks) == 1 and toks[0]["value"] == 1.2e-20
    flat = {"assoc.lead_p": 1.2e-20}
    v = verify_markdown("lead hit p = 1.2e-20", flat)
    assert v["passed"], v
    # a different exponent does NOT bind even at the same mantissa
    v2 = verify_markdown("lead hit p = 1.2e-8", flat)
    assert not v2["passed"]


def test_display_rounding_within_tolerance_binds():
    # recorded 0.25034 displayed as "25.0%" must still bind
    flat = {"frac": 0.25034}
    v = verify_markdown("25.0% of cells", flat)
    assert v["passed"], v


def test_binding_provenance_recorded():
    summary = _summary()
    md = render_markdown(summary)
    v = verify_markdown(md, _flat(summary), labels={"0", "1"})
    by_token = {c["token"]: c["bound_to"] for c in v["claims"]}
    # a waterfall number traces to the waterfall artifact...
    assert by_token["900"] == "cells_in"
    # ...and every bound claim points at a leaf that actually equals it
    flat = _flat(summary)
    for c in v["claims"]:
        if c["bound_to"] == "identifier":
            continue
        assert c["bound_to"] in flat, c
        val = float(c["token"].rstrip("%").replace(",", ""))
        if c["token"].endswith("%"):
            val /= 100.0
        assert abs(flat[c["bound_to"]] - val) < 0.51, c
    # a % cell binds the recorded fraction leaf
    assert by_token["44.96%"] == "cluster_qc.0.frac_of_total"


def test_cluster_labels_are_not_claims():
    flat = _flat(_summary())
    # '7' is a cluster label with no numeric leaf equal to 7
    v = verify_markdown("cluster 7 is the smallest", flat, labels={"7"})
    assert v["passed"], v
    assert v["claims"][0]["bound_to"] == "identifier"


def test_forbidden_phrase_flagged():
    v = verify_markdown("this proves the clusters are real",
                        _flat(_summary()))
    assert not v["passed"]
    assert "proves" in v["forbidden"]


# --- report() integration ---

def test_report_writes_verified_markdown(tmp_path, monkeypatch):
    import anndata as ad
    import numpy as np
    import scanpy as sc

    pytest.importorskip("scanpy")
    rng = np.random.default_rng(0)
    obs = pd.DataFrame({
        "leiden": pd.Categorical(["0"] * 60 + ["1"] * 40),
        "total_counts": rng.integers(500, 3000, 100).astype(float),
        "n_genes_by_counts": rng.integers(200, 900, 100).astype(float),
        "pct_mt": rng.uniform(0, 8, 100),
    })
    adata = ad.AnnData(X=rng.poisson(1, (100, 50)).astype(float), obs=obs)
    adata.obsm["X_umap"] = rng.normal(size=(100, 2))
    monkeypatch.setattr(sc.pl, "umap", lambda *a, **k: None)

    markers = pd.DataFrame({"group": ["0", "1"], "names": ["g1", "g2"],
                            "scores": [5.0, 4.0]})
    waterfall = {"cells_in": 120, "genes_in": 50,
                 "steps": [{"stage": "min_genes_per_cell",
                            "cells_before": 120, "cells_after": 100,
                            "cells_removed": 20, "genes_before": 50,
                            "genes_after": 50, "genes_removed": 0}],
                 "cells_out": 100, "genes_out": 50}
    cfg = {"dataset": {"mode": "demo"}}

    out = {k: str(tmp_path / k) for k in
           ("c.csv", "s.json", "u.png", "r.md", "claims.json")}
    monkeypatch.setattr("matplotlib.pyplot.savefig", lambda *a, **k: None)
    summary = report(adata, markers, waterfall, cfg,
                     out["c.csv"], out["s.json"], out["u.png"],
                     out["r.md"], out["claims.json"])

    md = (tmp_path / "r.md").read_text()
    assert "| cluster |" in md and f"{summary['n_cells']}" in md
    verdict = json.loads((tmp_path / "claims.json").read_text())
    assert verdict["passed"], verdict
    # mutating a number breaks the binding — the check must catch it
    tampered = md + "\n- total reads aligned: 424242\n"
    flat = {**flatten_results(waterfall), **flatten_results(summary)}
    labels = {r["cluster"] for r in summary["cluster_qc"]}
    assert not verify_markdown(tampered, flat, labels)["passed"]
