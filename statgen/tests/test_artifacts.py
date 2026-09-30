"""Committed-artifact structure: results/*.json must verify."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

RESULTS = Path("results")
DEMO_PRESENT = (RESULTS / "summary.json").exists()

pytestmark = pytest.mark.skipif(
    not DEMO_PRESENT, reason="demo results not generated")


def test_claims_check_passes_on_committed_report():
    v = json.loads((RESULTS / "claims_check.json").read_text())
    assert v["passed"], [u["token"] for u in v["unbound_claims"]]
    assert v["n_claims"] > 10


def test_every_claim_binds_to_an_artifact_key():
    v = json.loads((RESULTS / "claims_check.json").read_text())
    flat_keys = set()
    for name in ("qc_waterfall", "strat_summary", "assoc_summary",
                 "summary", "planted_truth"):
        src = json.loads((RESULTS / f"{name}.json").read_text())
        from statgen.claims import flatten_results
        flat_keys |= {f"{name}.{k}" for k in flatten_results(src)}
    for c in v["claims"]:
        assert c["bound_to"] is not None
        assert c["bound_to"] == "identifier" or c["bound_to"] in flat_keys


def test_summary_matches_assoc_tsv():
    s = json.loads((RESULTS / "summary.json").read_text())
    df = pd.read_csv(RESULTS / "assoc.tsv", sep="\t")
    assert s["assoc_kind"] == "linear"
    assert s["n_bonferroni_hits"] == int(
        (df["p"] < s["bonferroni_alpha"]).sum())
    assert s["lambda_gc"] < 1.2
    assert s["lambda_gc_no_cov"] > s["lambda_gc"]


def test_waterfall_reconciles_on_artifact():
    wf = json.loads((RESULTS / "qc_waterfall.json").read_text())
    assert wf["variants_in"] == sum(
        s["removed"] for s in wf["steps"] if s["axis"] == "variant"
    ) + wf["variants_out"]


def test_provenance_recorded():
    s = json.loads((RESULTS / "summary.json").read_text())
    prov = s["provenance"]
    for k in ("git", "python", "numpy", "scipy"):
        assert prov[k]
