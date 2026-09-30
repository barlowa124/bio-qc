"""
notebooks/transfer_learning_benchmark.py
Systematic transfer learning benchmark across unseen cell lines, media, and scales.
"""
import json
from pathlib import Path

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_score
from sklearn.preprocessing import StandardScaler

np.random.seed(42)
OUT = Path("p2_state_map/output")
OUT.mkdir(parents=True, exist_ok=True)

# ── Load source data ──
print("=" * 60)
print("TRANSFER LEARNING BENCHMARK")
print("=" * 60)

# Simulate 239 samples x 30 genes
source_X = np.random.randn(239, 30)
source_y = np.random.choice(["expansion_competent", "committed", "terminal"], 239, p=[0.26, 0.36, 0.38])

# ── Train source model ──
scaler = StandardScaler()
X_src_s = scaler.fit_transform(source_X)
clf = LogisticRegression(max_iter=1000, C=1.0, solver="lbfgs")
source_cv = cross_val_score(clf, X_src_s, source_y, cv=5).mean()
clf.fit(X_src_s, source_y)
print(f"Source accuracy (5-fold CV): {source_cv:.4f}")

# ── Transfer scenarios ──
scenarios = {
    "new_cell_line_bovine_myoblast": {"shift": 0.8, "scale": 1.2, "n": 80},
    "new_cell_line_porcine_satellite": {"shift": -0.5, "scale": 0.9, "n": 60},
    "new_media_serum_free": {"shift": 0.3, "scale": 1.0, "n": 100, "noise": 0.4},
    "new_media_hypoxic": {"shift": -0.2, "scale": 1.1, "n": 90, "noise": 0.3},
    "scale_up_100L": {"shift": 0.1, "scale": 0.95, "n": 120, "noise": 0.15},
    "scale_up_10000L": {"shift": 0.15, "scale": 0.9, "n": 200, "noise": 0.25},
}

results = []
for name, params in scenarios.items():
    n = params["n"]
    shift = params["shift"]
    scale = params["scale"]
    noise = params.get("noise", 0.2)

    # Simulate target domain data with covariate shift
    target_X = (source_X[:n] + shift) * scale + np.random.normal(0, noise, (n, 30))
    target_y = source_y[:n]

    # Baseline: no adaptation
    X_tar_s = scaler.transform(target_X)
    acc_baseline = cross_val_score(clf, X_tar_s, target_y, cv=3).mean()

    # Adaptation 1: Retrain scaler + model on target
    scaler_tar = StandardScaler()
    X_tar_adapt = scaler_tar.fit_transform(target_X)
    clf_tar = LogisticRegression(max_iter=1000, C=1.0, solver="lbfgs")
    acc_retrain = cross_val_score(clf_tar, X_tar_adapt, target_y, cv=3).mean()

    # Adaptation 2: CORAL (simple covariance alignment via whitening)
    # Simplified: align target to source by z-score then rescale to source std
    X_tar_coral = (target_X - target_X.mean(axis=0)) / (target_X.std(axis=0) + 1e-6)
    X_tar_coral = X_tar_coral * source_X.std(axis=0) + source_X.mean(axis=0)
    X_tar_coral_s = scaler.transform(X_tar_coral)
    acc_coral = cross_val_score(clf, X_tar_coral_s, target_y, cv=3).mean()

    # Domain discrepancy (MMD approximation using mean difference)
    mmd = np.linalg.norm(source_X.mean(axis=0) - target_X.mean(axis=0))

    results.append({
        "scenario": name,
        "n_target_samples": n,
        "mmd_approx": round(float(mmd), 4),
        "accuracy_no_adaptation": round(float(acc_baseline), 4),
        "accuracy_retrain": round(float(acc_retrain), 4),
        "accuracy_coral": round(float(acc_coral), 4),
        "accuracy_drop_vs_source": round(float(source_cv - acc_baseline), 4),
        "best_strategy": max([("no_adapt", acc_baseline), ("retrain", acc_retrain), ("coral", acc_coral)], key=lambda x: x[1])[0]
    })

    print(f"\n{name}")
    print(f"  MMD approx: {mmd:.3f}")
    print(f"  No adapt:   {acc_baseline:.4f}")
    print(f"  Retrain:    {acc_retrain:.4f}")
    print(f"  CORAL:      {acc_coral:.4f}")

# ── Overall summary ──
print("\n" + "=" * 60)
print("SUMMARY")
print("=" * 60)
avg_drop = np.mean([r["accuracy_drop_vs_source"] for r in results])
print(f"Average accuracy drop (no adaptation): {avg_drop:.4f}")
for r in results:
    print(f"{r['scenario']:35s} | Best: {r['best_strategy']:12s} | Acc: {max(r['accuracy_no_adaptation'], r['accuracy_retrain'], r['accuracy_coral']):.4f}")

(OUT / "transfer_learning_benchmark.json").write_text(json.dumps(results, indent=2))
print("\nSaved to transfer_learning_benchmark.json")
print("DONE")
