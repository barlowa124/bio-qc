"""Mondrian conformal on cluster-assignment confidence."""

import json

import numpy as np
import pytest

from scrna_qc.cluster_confidence import conformal_cluster_confidence


def _separable(n=800, k=3, d=10, seed=0):
    """Well-separated clusters in X — confidence should be high and
    coverage should meet the 1-alpha floor."""
    rng = np.random.default_rng(seed)
    labels = np.repeat(np.arange(k).astype(str), n // k)
    centers = rng.normal(size=(k, d)) * 5
    X = rng.normal(size=(len(labels), d)) + centers[labels.astype(int)]
    return X, labels


def test_pooled_coverage_meets_target_on_separable_data():
    pytest.importorskip("sklearn")
    X, labels = _separable()
    res = conformal_cluster_confidence(X, labels, alpha=0.1,
                                       cal_frac=0.3, min_group=30,
                                       seed=0)
    # finite-sample conformal: held-out coverage >= 1 - alpha, allowing
    # binomial slack on ~560 test cells
    assert res["pooled_coverage"] >= 0.85
    assert res["pooled_coverage"] <= 1.0
    assert res["n_cal"] + res["n_test"] == len(labels)


def test_per_cluster_coverage_reported_not_pooled_only():
    pytest.importorskip("sklearn")
    X, labels = _separable()
    res = conformal_cluster_confidence(X, labels, alpha=0.1,
                                       cal_frac=0.3, min_group=30,
                                       seed=0)
    assert set(res["per_cluster"]) == {"0", "1", "2"}
    for cl, row in res["per_cluster"].items():
        assert row["n_test"] > 0
        assert row["coverage"] is not None
        assert 0.0 <= row["coverage"] <= 1.0
        # sets can be empty on ambiguous cells — empty means the
        # evidence doesn't support ANY cluster at this level
        assert 0.0 <= row["mean_set_size"] <= 3.0


def test_small_cluster_falls_back_to_global_qhat():
    pytest.importorskip("sklearn")
    rng = np.random.default_rng(1)
    # 3 big clusters + one tiny cluster under min_group
    X_big, labels_big = _separable(n=750, k=3)
    X_small = rng.normal(size=(8, 10)) + 30
    X = np.vstack([X_big, X_small])
    labels = np.concatenate([labels_big, np.array(["tiny"] * 8)])
    res = conformal_cluster_confidence(X, labels, alpha=0.1,
                                       cal_frac=0.5, min_group=30,
                                       seed=1)
    assert res["per_cluster"]["tiny"]["global_fallback"] is True
    assert res["per_cluster"]["tiny"]["qhat"] == res["qhat_global"]
    for cl in ("0", "1", "2"):
        assert res["per_cluster"][cl]["global_fallback"] is False


def test_confident_clusters_have_small_sets():
    pytest.importorskip("sklearn")
    X, labels = _separable()
    res = conformal_cluster_confidence(X, labels, alpha=0.1,
                                       cal_frac=0.3, min_group=30,
                                       seed=0)
    # separable data: near-singleton sets are the honest outcome
    assert res["mean_set_size"] < 1.5


def test_deterministic_given_seed():
    pytest.importorskip("sklearn")
    X, labels = _separable(n=300, k=3)
    a = conformal_cluster_confidence(X, labels, 0.1, 0.3, 30, 0)
    b = conformal_cluster_confidence(X, labels, 0.1, 0.3, 30, 0)
    assert a == b
