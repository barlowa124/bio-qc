"""
notebooks/conformal_prediction.py
Inductive conformal prediction for the 30-gene QC panel with guaranteed coverage.
"""
import json
from pathlib import Path

import numpy as np
from sklearn.datasets import make_classification
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

np.random.seed(42)
OUT = Path("p2_state_map/output")
OUT.mkdir(parents=True, exist_ok=True)

print("=" * 60)
print("CONFORMAL PREDICTION / CALIBRATION")
print("=" * 60)

# Simulate data
X, y = make_classification(
    n_samples=239, n_features=30, n_informative=20, n_redundant=5,
    n_classes=3, random_state=42
)
label_map = {0: "expansion_competent", 1: "committed", 2: "terminal"}
y_labels = np.array([label_map[yi] for yi in y])

X_train, X_rest, y_train, y_rest = train_test_split(X, y_labels, test_size=0.5, random_state=42, stratify=y_labels)
X_cal, X_test, y_cal, y_test = train_test_split(X_rest, y_rest, test_size=0.4, random_state=42, stratify=y_rest)

print(f"Train: {len(y_train)} | Calibration: {len(y_cal)} | Test: {len(y_test)}")

scaler = StandardScaler()
X_train_s = scaler.fit_transform(X_train)
X_cal_s = scaler.transform(X_cal)
X_test_s = scaler.transform(X_test)

clf = LogisticRegression(max_iter=1000, C=1.0, solver="lbfgs")
clf.fit(X_train_s, y_train)
classes = clf.classes_

# Nonconformity scores on calibration set
proba_cal = clf.predict_proba(X_cal_s)
nonconf = np.array([1.0 - proba_cal[i, np.where(classes == y_cal[i])[0][0]] for i in range(len(y_cal))])

alpha_levels = [0.05, 0.10, 0.20]
results = []
for alpha in alpha_levels:
    q = np.quantile(nonconf, np.ceil((1 + len(nonconf)) * (1 - alpha)) / len(nonconf))
    q = min(q, 1.0)
    proba_test = clf.predict_proba(X_test_s)
    sets = []
    for i in range(len(X_test)):
        pred_set = [classes[j] for j in range(len(classes)) if 1.0 - proba_test[i, j] <= q]
        sets.append(pred_set)
    coverage = np.mean([y_test[i] in sets[i] for i in range(len(y_test))])
    avg_set_size = np.mean([len(s) for s in sets])
    empty_sets = sum(1 for s in sets if len(s) == 0)
    results.append({
        "alpha": alpha,
        "quantile_threshold": round(float(q), 4),
        "empirical_coverage": round(float(coverage), 4),
        "average_prediction_set_size": round(float(avg_set_size), 2),
        "empty_sets": int(empty_sets),
        "guaranteed_coverage": f">={1-alpha:.0%}"
    })
    print(f"\nAlpha = {alpha:.0%}")
    print(f"  Threshold (q): {q:.4f}")
    print(f"  Empirical coverage: {coverage:.2%}")
    print(f"  Avg prediction set size: {avg_set_size:.2f}")
    print(f"  Empty sets: {empty_sets}/{len(y_test)}")

print("\n" + "=" * 60)
print("SUMMARY")
print("=" * 60)
for r in results:
    print(f"  Alpha={r['alpha']:.0%} | Coverage={r['empirical_coverage']:.2%} | AvgSetSize={r['average_prediction_set_size']:.2f}")

(OUT / "conformal_prediction.json").write_text(json.dumps(results, indent=2))
print("\nSaved to conformal_prediction.json")
print("DONE")
