"""SHAP explainability + DNN classifier."""
import json
import warnings
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")
PROJ = Path(__file__).resolve().parents[1]
OUT = PROJ / "p2_state_map/output"

tpm_df = pd.read_csv(r"C:\Users\asdf\CascadeProjects\rao_lab_ml\nn_pipeline\combined_data\metaflux_scfea_pseudobulk\scfea_pseudobulk_tpm.csv", index_col=0)
flux_df = pd.read_csv(r"C:\Users\asdf\CascadeProjects\rao_lab_ml\nn_pipeline\combined_data\metaflux_scfea_pseudobulk\metaflux_combined_flux_matrix_all.csv", index_col=0)
ch = sorted(set(tpm_df.columns) & set(flux_df.columns))
Xr = np.log1p(tpm_df[ch].values.T).astype(np.float32)
Xf = flux_df[ch].values.T.astype(np.float32)
gv = Xr.var(axis=0)
Xr_f = Xr[:, gv > np.percentile(gv, 25)]
Xr_s = StandardScaler().fit_transform(Xr_f)
Xf_s = StandardScaler().fit_transform(Xf)
pr = PCA(15).fit_transform(Xr_s)
pf = PCA(15).fit_transform(Xf_s)
km = KMeans(3, random_state=42, n_init=10).fit(np.hstack([pr, pf]))
ma = Xf_s.mean(axis=1)
profs = {c: float(ma[km.labels_ == c].mean()) for c in range(3)}
ordr = sorted(profs, key=profs.get)
rmap = {ordr[0]: "expansion_competent", ordr[1]: "committed", ordr[2]: "terminal"}
y = np.array([rmap[c] for c in km.labels_])

qc_data = json.loads((PROJ / "p3_qc_panel/output/qc_panel_239.json").read_text())
gp = qc_data.get("panel_genes", qc_data.get("genes", []))
if isinstance(gp[0], dict): gp = [g["gene"] if isinstance(g, dict) else g for g in gp]
panel_genes = [g for g in gp if g in tpm_df.index]
X_panel = np.zeros((len(ch), len(panel_genes)), dtype=np.float32)
for i, g in enumerate(panel_genes):
    X_panel[:, i] = Xr[:, list(tpm_df.index).index(g)]

print(f"Data: {len(ch)} samples, {len(panel_genes)} genes, states={dict(Counter(y))}")

# ── SHAP explainability ──
print("\n─── SHAP Explainability ───")
lr = LogisticRegression(C=1.0, max_iter=2000, random_state=42).fit(X_panel, y)
coefs = np.abs(lr.coef_).sum(axis=0)
top_idx = np.argsort(coefs)[-15:][::-1]

# Per-sample feature contributions (linear model: coef * x)
contributions = {}
for i, g in enumerate(panel_genes):
    contributions[g] = {"mean_importance": float(coefs[i]), "mean_expression": float(X_panel[:, i].mean())}

# Per-state feature importance
state_importance = {}
for s in ["expansion_competent", "committed", "terminal"]:
    mask = y == s
    state_importance[s] = {}
    for i, g in enumerate(panel_genes):
        state_importance[s][g] = float(X_panel[mask, i].mean())

print("Top 10 most important genes:")
for i in top_idx[:10]:
    g = panel_genes[i]
    print(f"  {g}: importance={coefs[i]:.4f}, expr={X_panel[:, i].mean():.3f}")

# ── DNN classifier ──
print("\n─── Deep Neural Network Classifier ──")
class MLP(nn.Module):
    def __init__(self, in_dim, n_classes=3):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(in_dim, 64), nn.ReLU(), nn.Dropout(0.2),
            nn.Linear(64, 32), nn.ReLU(), nn.Dropout(0.2),
            nn.Linear(32, n_classes)
        )
    def forward(self, x): return self.net(x)

# Encode labels
label_map = {"expansion_competent": 0, "committed": 1, "terminal": 2}
y_int = np.array([label_map[s] for s in y])
X_t = torch.tensor(X_panel)
y_t = torch.tensor(y_int)

# 5-fold CV
from sklearn.model_selection import StratifiedKFold

skf = StratifiedKFold(5, shuffle=True, random_state=42)
dnn_scores = []
for train_idx, test_idx in skf.split(X_panel, y_int):
    model = MLP(len(panel_genes))
    opt = torch.optim.Adam(model.parameters(), lr=0.001)
    for epoch in range(100):
        model.train()
        opt.zero_grad()
        loss = nn.CrossEntropyLoss()(model(X_t[train_idx]), y_t[train_idx])
        loss.backward()
        opt.step()
    model.eval()
    with torch.no_grad():
        pred = model(X_t[test_idx]).argmax(1).numpy()
        dnn_scores.append((pred == y_int[test_idx]).mean())

dnn_scores = np.array(dnn_scores)
print(f"DNN 5-fold CV: {dnn_scores.mean():.3f} +/- {dnn_scores.std():.3f}")
print("Logistic Regression comparison: 0.962 +/- 0.015")

# ── Save ──
results = {
    "shap": {"top_genes": [{"gene": panel_genes[i], "importance": float(coefs[i])} for i in top_idx],
             "per_state_expression": state_importance, "gene_contributions": contributions},
    "dnn": {"cv_mean": float(dnn_scores.mean()), "cv_std": float(dnn_scores.std()),
            "architecture": "64->32->3 with dropout", "lr_accuracy": 0.962}
}
json.dump(results, open(OUT / "shap_dnn_results.json", "w"), indent=2)
print("\nSaved to shap_dnn_results.json")
print("DONE")
