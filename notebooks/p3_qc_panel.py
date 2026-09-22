"""P3: Minimal QC Biomarker Panel for Early Batch Triage."""
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from sklearn.feature_selection import SelectFromModel
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import LeaveOneGroupOut, cross_val_score
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")
PROJ = Path(__file__).resolve().parents[1]
CFG = yaml.safe_load((PROJ / "shared_data/data_paths.yaml").read_text())
OUT = PROJ / "p3_qc_panel/output"
OUT.mkdir(exist_ok=True)

print("=" * 60)
print("P3: MINIMAL QC BIOMARKER PANEL")
print("=" * 60)

# 1. Load
print("\n1. Loading data...")
norm = pd.read_csv(CFG["rna_seq"]["normalized_counts"], index_col=0)
de = pd.read_csv(CFG["rna_seq"]["de_results"])
meta = pd.read_csv(CFG["metadata"]["geo"])
conditions = meta["Condition"].values

X = np.log1p(norm.values.T)
gene_names = norm.index.values
print(f"   Expression: {X.shape}, DE genes: {len(de)}")

# 2. Define QC labels
print("\n2. Defining QC labels...")
# TSCM = good (expansion-competent stem cells)
# TA = acceptable (committed but viable)
# TUA = needs investigation (terminal, may indicate drift)
label_map = {}
for c in meta["Condition"].unique():
    if "TSCM" in c:
        label_map[c] = "pass"
    elif "TA" in c:
        label_map[c] = "pass"
    elif "TUA" in c:
        label_map[c] = "review"
qc_labels = np.array([label_map[c] for c in conditions])
print(f"   Pass: {(qc_labels == 'pass').sum()}, Review: {(qc_labels == 'review').sum()}")

# 3. Feature selection: top DE genes
print("\n3. Selecting candidate biomarkers...")
# Use DE results to narrow to significant genes
de_sig = de[de["padj"].notna() & (de["padj"] < 0.05)]
top_de = de_sig.nlargest(200, "abs(log2FoldChange)") if "abs(log2FoldChange)" in de_sig.columns else de_sig.nlargest(200, "stat")
candidate_genes = [g for g in top_de["Gene"].values if g in gene_names]
print(f"   Candidate genes (sig DE + in expression matrix): {len(candidate_genes)}")

# Build candidate matrix
cand_idx = [list(gene_names).index(g) for g in candidate_genes]
X_cand = X[:, cand_idx]
X_cand_s = StandardScaler().fit_transform(X_cand)

# 4. L1-regularized logistic regression for sparse selection
print("\n4. L1 feature selection...")
groups = np.array([c.split("_")[0] for c in conditions])
logo = LeaveOneGroupOut()

lr_l1 = LogisticRegression(penalty="l1", solver="saga", max_iter=5000, C=0.1)
selector = SelectFromModel(lr_l1, max_features=30, threshold=-np.inf)
selector.fit(X_cand_s, qc_labels)
selected_mask = selector.get_support()
selected_genes = [candidate_genes[i] for i in range(len(candidate_genes)) if selected_mask[i]]
print(f"   Selected genes: {len(selected_genes)}")
for g in selected_genes:
    print(f"     {g}")

# 5. Cross-validate
print("\n5. Cross-validation...")
X_sel = X_cand_s[:, selected_mask]
lr_final = LogisticRegression(max_iter=1000)
scores = cross_val_score(lr_final, X_sel, qc_labels, cv=logo, groups=groups)
print(f"   Leave-one-sex-out accuracy: {scores.mean():.2f} (+/- {scores.std():.2f})")

# 6. Decision rule
print("\n6. Decision rule...")
lr_final.fit(X_sel, qc_labels)
print(f"   Classes: {lr_final.classes_.tolist()}")
print(f"   Intercept: {lr_final.intercept_.tolist()}")
for i, g in enumerate(selected_genes[:10]):
    print(f"   {g}: coef={lr_final.coef_[0][i]:.4f}")

# 7. Save
panel = {
    "n_genes": len(selected_genes),
    "genes": selected_genes,
    "accuracy": float(scores.mean()),
    "accuracy_std": float(scores.std()),
    "decision_classes": lr_final.classes_.tolist(),
    "intercept": lr_final.intercept_.tolist(),
    "coefficients": {g: float(c) for g, c in zip(selected_genes, lr_final.coef_[0])},
}
with open(OUT / "qc_panel.json", "w") as f:
    json.dump(panel, f, indent=2)

print(f"\nOutputs saved to: {OUT}")
print("=" * 60)
print("P3 COMPLETE")
print("=" * 60)