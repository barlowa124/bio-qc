"""Cluster-assignment confidence with Mondrian conformal coverage.

Leiden labels are hard assignments; this stage gives them a calibrated
set-valued form. A logistic classifier on the PCA coordinates produces
per-cluster probabilities; the shared split-conformal helpers (vendored
conformal.py — fourth consumer after comp-tox, protstab, and
cultivated-meat) calibrate nonconformity per true cluster, falling back
to the pooled level for clusters below `min_group` calibration cells.

The reported guarantee is conditional, not pooled: a pooled 90%
prediction-set coverage can sit at 60% inside a rare cluster, which is
exactly the failure mode per-cluster QC exists to expose. Coverage
numbers on the held-out split are empirical measurements, not the
theoretical bound — small clusters are noisier, and `min_group` marks
which clusters used the global fallback.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np

from .conformal import class_scores, covered_set, mondrian_qhat
from .util import load_config


def conformal_cluster_confidence(X: np.ndarray, labels: np.ndarray,
                                 alpha: float, cal_frac: float,
                                 min_group: int, seed: int) -> dict:
    """Split cal/test, fit LR on X, mondrian-calibrate, measure coverage."""
    from sklearn.linear_model import LogisticRegression

    rng = np.random.default_rng(seed)
    classes = np.unique(labels)
    cal_idx, test_idx = [], []
    for c in classes:
        idx = np.flatnonzero(labels == c)
        rng.shuffle(idx)
        k = max(1, int(round(len(idx) * cal_frac)))
        cal_idx += idx[:k].tolist()
        test_idx += idx[k:].tolist()
    cal_idx, test_idx = np.array(cal_idx), np.array(test_idx)

    clf = LogisticRegression(max_iter=1000, random_state=seed)
    clf.fit(X[cal_idx], labels[cal_idx])
    class_to_col = {c: i for i, c in enumerate(clf.classes_)}
    y_cal = np.array([class_to_col[l] for l in labels[cal_idx]])
    y_test = np.array([class_to_col[l] for l in labels[test_idx]])

    probs_cal = clf.predict_proba(X[cal_idx])
    probs_test = clf.predict_proba(X[test_idx])
    scores = class_scores(probs_cal, y_cal)
    qs = mondrian_qhat(scores, labels[cal_idx], alpha, min_group)

    # mondrian prediction set: include class k when 1 - p_k <= q[k]
    sets = np.zeros(probs_test.shape, dtype=bool)
    for c, col in class_to_col.items():
        sets[:, col] = probs_test[:, col] >= (1.0 - qs[c])

    covered = covered_set(y_test, sets)
    per_cluster = {}
    for c in classes:
        sel = labels[test_idx] == c
        n = int(sel.sum())
        per_cluster[str(c)] = {
            "n_test": n,
            "coverage": (float(covered[sel].mean()) if n else None),
            "qhat": qs[c],
            "global_fallback": len(np.flatnonzero(labels[cal_idx] == c))
                               < min_group,
            "mean_set_size": (float(sets[sel].sum(1).mean()) if n else None),
        }

    return {
        "alpha": alpha,
        "n_cal": int(len(cal_idx)),
        "n_test": int(len(test_idx)),
        "pooled_coverage": float(covered.mean()),
        "mean_set_size": float(sets.sum(1).mean()),
        "qhat_global": qs["global"],
        "min_group": min_group,
        "per_cluster": per_cluster,
        "estimator": "LogisticRegression(max_iter=1000) on X_pca",
    }


def main() -> None:
    emb_in, out_json = sys.argv[1:3]
    cfg = load_config()
    cc = cfg["conformal"]
    adata = ad.read_h5ad(emb_in)
    X = np.asarray(adata.obsm["X_pca"])
    labels = np.asarray(adata.obs["leiden"].astype(str))
    res = conformal_cluster_confidence(
        X, labels, alpha=cc["alpha"], cal_frac=cc["cal_frac"],
        min_group=cc["min_group"], seed=int(cfg["embed"]["seed"]))
    res["provenance"] = {
        "classifier": "sklearn.LogisticRegression",
        "features": "X_pca (post-HVG PCA)",
    }
    Path(out_json).parent.mkdir(parents=True, exist_ok=True)
    with open(out_json, "w") as f:
        json.dump(res, f, indent=2)
    print(f"conformal: pooled coverage {res['pooled_coverage']:.3f} at "
          f"alpha={res['alpha']} across {len(res['per_cluster'])} clusters")


if __name__ == "__main__":
    main()
