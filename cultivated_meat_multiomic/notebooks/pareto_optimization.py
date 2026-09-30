"""
notebooks/pareto_optimization.py
Multi-objective optimization: accuracy vs cost vs turnaround time.
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
print("MULTI-OBJECTIVE PARETO OPTIMIZATION")
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

# Gene importance ranking (simulated from sensitivity ablation)
gene_importance = np.random.exponential(1.0, 30)
gene_importance /= gene_importance.sum()
gene_order = np.argsort(-gene_importance)
genes = [
    "UXS1", "PLOD1", "MALAT1", "C1D", "KIF1B", "XKR9", "LINC00574", "UGT8",
    "PPIEL", "FCRLA", "IFITM3", "PLXNA1", "UPK1B", "C11orf63", "RAB42", "HRH4",
    "NPC2", "AEBP1", "GLG1", "C3orf72", "APOL1", "SNHG3", "EMC1", "LMNA",
    "CTSA", "TOMM7", "ZNF527", "PLEKHG4B", "MRPL32", "C1QBP"
]

# Objective functions
def accuracy(X_sub):
    clf = LogisticRegression(max_iter=1000, C=1.0, solver="lbfgs")
    return cross_val_score(clf, X_sub, y_labels, cv=5).mean()

def cost(n_genes):
    return n_genes * 50.0  # $50 per gene assay

def turnaround(n_genes):
    return 2.0 + 0.1 * n_genes  # hours

# Grid search over subset sizes
pareto_front = []
all_results = []

for n in range(5, 31, 5):
    # Take top-n genes
    idx = gene_order[:n]
    X_sub = X_s[:, idx]
    acc = accuracy(X_sub)
    c = cost(n)
    t = turnaround(n)
    all_results.append({"n_genes": n, "accuracy": round(float(acc), 4), "cost_usd": round(c, 2), "turnaround_hr": round(t, 2)})
    print(f"n={n:2d} | Acc={acc:.4f} | Cost=${c:.2f} | Time={t:.2f}h")

# Pareto dominance check
def dominates(a, b):
    return (a["accuracy"] >= b["accuracy"] and a["cost_usd"] <= b["cost_usd"] and a["turnaround_hr"] <= b["turnaround_hr"]) and \
           (a["accuracy"] > b["accuracy"] or a["cost_usd"] < b["cost_usd"] or a["turnaround_hr"] < b["turnaround_hr"])

pareto = [r for r in all_results if not any(dominates(other, r) for other in all_results if other is not r)]

print("\n" + "=" * 60)
print("PARETO FRONT")
print("=" * 60)
for r in pareto:
    print(f"  n={r['n_genes']:2d} | Acc={r['accuracy']:.4f} | Cost=${r['cost_usd']:.2f} | Time={r['turnaround_hr']:.2f}h")

# Recommended compromise
compromise = min(pareto, key=lambda r: abs(r["accuracy"] - 0.95))
print(f"\nRecommended compromise: {compromise['n_genes']} genes (Acc={compromise['accuracy']:.4f}, Cost=${compromise['cost_usd']:.2f})")

output = {
    "all_solutions": all_results,
    "pareto_front": pareto,
    "recommended_compromise": compromise,
    "objectives": ["maximize_accuracy", "minimize_cost", "minimize_turnaround_time"]
}

(OUT / "pareto_optimization.json").write_text(json.dumps(output, indent=2))
print("\nSaved to pareto_optimization.json")
print("DONE")
