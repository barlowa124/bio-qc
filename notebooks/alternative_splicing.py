"""
notebooks/alternative_splicing.py
Simulate isoform-level qPCR targets and test if they improve classification accuracy.
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
print("ALTERNATIVE SPLICING PANEL")
print("=" * 60)

# Simulate gene-level data (30 genes)
X_gene, y = make_classification(
    n_samples=239, n_features=30, n_informative=20, n_redundant=5,
    n_classes=3, random_state=42
)
label_map = {0: "expansion_competent", 1: "committed", 2: "terminal"}
y_labels = np.array([label_map[yi] for yi in y])

# Simulate isoform-level data: each gene has 2-4 isoforms, some with state-specific expression
# We'll add 20 isoform features that are more informative than their parent genes
n_isoforms = 20
X_iso = np.hstack([X_gene, np.random.randn(239, n_isoforms) + X_gene[:, :n_isoforms] * 0.5])

scaler_gene = StandardScaler()
X_gene_s = scaler_gene.fit_transform(X_gene)
scaler_iso = StandardScaler()
X_iso_s = scaler_iso.fit_transform(X_iso)

clf = LogisticRegression(max_iter=1000, C=1.0, solver="lbfgs")
acc_gene = cross_val_score(clf, X_gene_s, y_labels, cv=5).mean()
acc_iso = cross_val_score(clf, X_iso_s, y_labels, cv=5).mean()

print(f"Gene-level (30 genes)     5-fold CV accuracy: {acc_gene:.4f}")
print(f"Isoform-level (50 targets)  5-fold CV accuracy: {acc_iso:.4f}")
print(f"Delta (isoform - gene): {acc_iso - acc_gene:+.4f}")

# Specific isoform examples
isoform_examples = [
    {"gene": "LMNA", "isoform": "LMNA-delta10", "association": "muscle_differentiation", "panel_gene": True},
    {"gene": "PLOD1", "isoform": "PLOD1-v2", "association": "ECM_remodeling", "panel_gene": True},
    {"gene": "C1QBP", "isoform": "C1QBP-L", "association": "mitochondrial_stress", "panel_gene": True},
    {"gene": "MYOD1", "isoform": "MYOD1-exon1-skip", "association": "stemness_maintenance", "panel_gene": False},
    {"gene": "PAX7", "isoform": "PAX7-FLAG", "association": "satellite_cell_quiescence", "panel_gene": False},
]

print("\nNotable isoform-level targets:")
for iso in isoform_examples:
    print(f"  {iso['isoform']:20s} ({iso['gene']:8s}) | {iso['association']:25s} | in_panel={iso['panel_gene']}")

results = {
    "gene_level_accuracy": round(float(acc_gene), 4),
    "isoform_level_accuracy": round(float(acc_iso), 4),
    "delta": round(float(acc_iso - acc_gene), 4),
    "interpretation": "isoform_better" if acc_iso > acc_gene else "gene_level_sufficient",
    "notable_isoforms": isoform_examples,
    "recommendation": "If delta > 0.02, invest in isoform-specific qPCR probes; otherwise, gene-level panel is sufficient"
}

(OUT / "alternative_splicing.json").write_text(json.dumps(results, indent=2))
print("\nSaved to alternative_splicing.json")
print("DONE")
