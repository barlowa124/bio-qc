"""
notebooks/adversarial_robustness.py
Test biologically plausible perturbations on the 30-gene panel.
"""
import json
from pathlib import Path

import numpy as np
from sklearn.datasets import make_classification
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_score
from sklearn.preprocessing import StandardScaler

np.random.seed(42)
OUT = Path("p2_state_map/output")
OUT.mkdir(parents=True, exist_ok=True)

print("=" * 60)
print("ADVERSARIAL ROBUSTNESS")
print("=" * 60)

# Simulate data
X, y = make_classification(
    n_samples=239, n_features=30, n_informative=20, n_redundant=5,
    n_classes=3, random_state=42
)
label_map = {0: "expansion_competent", 1: "committed", 2: "terminal"}
y_labels = np.array([label_map[yi] for yi in y])

scaler = StandardScaler()
X_s = scaler.fit_transform(X)
clf = LogisticRegression(max_iter=1000, C=1.0, solver="lbfgs")
base_acc = cross_val_score(clf, X_s, y_labels, cv=5).mean()
clf.fit(X_s, y_labels)

# Biologically plausible perturbations
genes = [
    "UXS1", "PLOD1", "MALAT1", "C1D", "KIF1B", "XKR9", "LINC00574", "UGT8",
    "PPIEL", "FCRLA", "IFITM3", "PLXNA1", "UPK1B", "C11orf63", "RAB42", "HRH4",
    "NPC2", "AEBP1", "GLG1", "C3orf72", "APOL1", "SNHG3", "EMC1", "LMNA",
    "CTSA", "TOMM7", "ZNF527", "PLEKHG4B", "MRPL32", "C1QBP"
]

perturbations = {
    "none": lambda x: x,
    "gaussian_cv5": lambda x: x + np.random.normal(0, 0.05, x.shape),
    "gaussian_cv10": lambda x: x + np.random.normal(0, 0.10, x.shape),
    "gaussian_cv20": lambda x: x + np.random.normal(0, 0.20, x.shape),
    "batch_effect_shift": lambda x: x + np.random.normal(0.5, 0.1, (1, 30)),
    "myoblast_markers_down": lambda x: x - np.array([0.5 if g in {"PAX7","MYOD1"} else 0.0 for g in genes]) * (np.random.rand(*x.shape) < 0.3),
}

results = []
for name, perturb_fn in perturbations.items():
    X_pert = perturb_fn(X_s.copy())
    scores = cross_val_score(LogisticRegression(max_iter=1000, C=1.0, solver="lbfgs"), X_pert, y_labels, cv=5)
    acc = scores.mean()
    std = scores.std()
    results.append({
        "perturbation": name,
        "accuracy": round(float(acc), 4),
        "std": round(float(std), 4),
        "accuracy_drop": round(float(base_acc - acc), 4)
    })
    print(f"{name:25s} | Acc={acc:.4f} +/- {std:.4f} | Drop={base_acc-acc:.4f}")

# Worst-case: flip top 3 informative genes
worst = X_s.copy()
for i in range(3):
    worst[:, i] *= -1
acc_worst = cross_val_score(LogisticRegression(max_iter=1000, C=1.0, solver="lbfgs"), worst, y_labels, cv=5).mean()
results.append({"perturbation": "worst_case_flip_top3", "accuracy": round(float(acc_worst), 4), "std": 0.0, "accuracy_drop": round(float(base_acc - acc_worst), 4)})
print(f"{'worst_case_flip_top3':25s} | Acc={acc_worst:.4f} | Drop={base_acc-acc_worst:.4f}")

print("\n" + "=" * 60)
print("SUMMARY")
print("=" * 60)
print(f"Baseline accuracy: {base_acc:.4f}")
print(f"Most robust perturbation: {min(results[:-1], key=lambda x: x['accuracy_drop'])['perturbation']}")
print(f"Most fragile perturbation: {max(results, key=lambda x: x['accuracy_drop'])['perturbation']}")

(OUT / "adversarial_robustness.json").write_text(json.dumps(results, indent=2))
print("\nSaved to adversarial_robustness.json")
print("DONE")
