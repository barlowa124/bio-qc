"""Unit tests on the math plus smoke tests on the generated cohort."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from statgen import assoc, qc, strat
from statgen.data import generate_demo

RESULTS = Path("results")


_DEMO_CACHE = None


def _demo():
    global _DEMO_CACHE
    if _DEMO_CACHE is None:
        cfg = {
            "dataset": {"demo": {
                "seed": 11, "n_samples": 800, "n_variants": 4000,
                "ancestry_shares": [0.6, 0.4], "fst": 0.08,
                "missing_rate": 0.02, "n_causal": 8,
                "beta_lo": 0.35, "beta_hi": 0.60,
                "ancestry_shift": 0.8, "h2": 0.35}}}
        _DEMO_CACHE = generate_demo(cfg)
    return _DEMO_CACHE


_QC_CACHE = None


def _qc():
    global _QC_CACHE
    if _QC_CACHE is None:
        cohort, _ = _demo()
        q = {"max_sample_missing": 0.10, "max_variant_missing": 0.05,
             "min_maf": 0.01, "hwe_p_min": 1e-6}
        _QC_CACHE = qc.apply_qc(cohort["G"], cohort["samples"], q)
    return _QC_CACHE


def test_demo_generator_is_deterministic():
    a, _ = _demo()
    b, _ = _demo()
    np.testing.assert_array_equal(a["G"], b["G"])
    np.testing.assert_allclose(a["y"], b["y"])


def test_planted_truth_records_effects():
    _, truth = _demo()
    assert truth["n_causal"] == len(truth["causal"]) == 8
    betas = [abs(c["beta"]) for c in truth["causal"]]
    assert all(0.35 <= b <= 0.60 for b in betas)


def test_maf_math():
    G = np.array([[0, 1, 2], [2, 2, 2]], dtype=np.int8)
    m = qc.maf(G)
    assert m.shape == (3,)
    np.testing.assert_allclose(m[2], 0.0)  # monomorphic
    np.testing.assert_allclose(m[0], 0.5)   # freq 0.5 -> maf 0.5
    np.testing.assert_allclose(m[1], 0.25)  # freq 0.75 -> maf 0.25


def test_hwe_detects_hard_departure():
    rng = np.random.default_rng(0)
    n = 500
    # all heterozygotes: maximal HWE violation at p=0.5
    G_bad = np.ones((n, 1), dtype=np.int8)
    assert qc.hwe_p(G_bad)[0] < 1e-20
    # HWE-distributed genotypes at p=0.5
    g = np.repeat([0, 1, 2], [125, 250, 125])
    rng.shuffle(g)
    assert qc.hwe_p(g.reshape(-1, 1))[0] > 0.9


def test_waterfall_reconciles():
    _, _, _, _, wf = _qc()
    assert wf["variants_in"] == sum(
        s["removed"] for s in wf["steps"] if s["axis"] == "variant"
    ) + wf["variants_out"]
    assert wf["samples_in"] == wf["samples_out"] + sum(
        s["removed"] for s in wf["steps"] if s["axis"] == "sample")


def test_pc1_recovers_ancestry_structure():
    cohort, _ = _demo()
    G, _, keep_s, _, _ = _qc()
    pcs, eig = strat.run(G, 4)
    anc = cohort["ancestry"][keep_s].astype(float)
    r = abs(np.corrcoef(pcs[:, 0], anc)[0, 1])
    assert r > 0.9, f"PC1 should track ancestry, got r={r:.2f}"


def _assoc_inputs():
    cohort, truth = _demo()
    G, _, keep_s, idx_v, _ = _qc()
    Gf = assoc._impute(G)
    y = cohort["y"][keep_s]
    pcs, _ = strat.run(G, 2)
    X = np.column_stack(
        [np.ones(len(y)), cohort["cov"][keep_s], pcs])
    return Gf, y, X, idx_v, truth, cohort


def test_linear_scan_recovers_planted_betas():
    Gf, y, X, idx_v, truth, cohort = _assoc_inputs()
    beta, _, stat, p = assoc.linear_scan(y, Gf, X)
    causal = {c["vid"]: c for c in truth["causal"]}
    vids = cohort["vid"][idx_v]
    hits = 0
    for i, v in enumerate(vids):
        if v in causal:
            if p[i] < 0.05 / len(vids):
                hits += 1
            assert beta[i] == pytest.approx(causal[v]["beta"], abs=0.2)
    assert hits >= 3


def test_pc_covariates_control_inflation():
    Gf, y, X, _, _, _ = _assoc_inputs()
    _, _, stat_cov, _ = assoc.linear_scan(y, Gf, X)
    _, _, stat_null, _ = assoc.linear_scan(
        y, Gf, np.ones((len(y), 1)))
    lam_cov = assoc.genomic_lambda(stat_cov)
    lam_null = assoc.genomic_lambda(stat_null)
    assert lam_null > 1.3, "stratified cohort should inflate"
    assert lam_cov < lam_null
    assert lam_cov < 1.2, "PC covariates should control inflation"


def test_logistic_score_scan_smoke():
    Gf, y, X, _, _, _ = _assoc_inputs()
    yb = (y > np.median(y)).astype(float)
    beta, se, stat, p = assoc.logistic_score_scan(yb, Gf, X)
    assert p.shape == (Gf.shape[1],)
    assert np.all((p >= 0) & (p <= 1))
    lam = assoc.genomic_lambda(stat)
    assert 0.5 < lam < 1.5


def test_missing_imputation_preserves_shape():
    G = np.array([[0, -1], [2, 1], [-1, 2]], dtype=np.int8)
    Gi = assoc._impute(G)
    assert Gi.shape == G.shape and not np.isnan(Gi).any()
    assert Gi[0, 1] == pytest.approx(1.5)
