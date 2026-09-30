"""Remaining analyses: SHAP, WGCNA, cell communication, multi-omic benchmark, DNN, DE, batch effects."""
import json
import warnings
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")
PROJ = Path(__file__).resolve().parents[1]
CROSS = PROJ / "cross_species_validation"
OUT = PROJ / "p2_state_map/output"
OUT.mkdir(exist_ok=True)

print("=" * 70)
print("REMAINING COMPUTATIONAL ANALYSES")
print("=" * 70)

# ═══════════════════════════════════════════════════════════
# SHARED DATA
# ═══════════════════════════════════════════════════════════

tpm_df = pd.read_csv(r"C:\Users\asdf\CascadeProjects\rao_lab_ml\nn_pipeline\combined_data\metaflux_scfea_pseudobulk\scfea_pseudobulk_tpm.csv", index_col=0)
flux_df = pd.read_csv(r"C:\Users\asdf\CascadeProjects\rao_lab_ml\nn_pipeline\combined_data\metaflux_scfea_pseudobulk\metaflux_combined_flux_matrix_all.csv", index_col=0)
ch = sorted(set(tpm_df.columns) & set(flux_df.columns))
Xr = np.log1p(tpm_df[ch].values.T).astype(np.float32)
Xf = flux_df[ch].values.T.astype(np.float32)
gv = Xr.var(axis=0)
kp = gv > np.percentile(gv, 25)
Xr_f = Xr[:, kp]
Xr_s = StandardScaler().fit_transform(Xr_f)
Xf_s = StandardScaler().fit_transform(Xf)
pr = PCA(n_components=15).fit_transform(Xr_s)
pf = PCA(n_components=15).fit_transform(Xf_s)
Xj = np.hstack([pr, pf])
km = KMeans(n_clusters=3, random_state=42, n_init=10)
cl = km.fit_predict(Xj)
ma = Xf_s.mean(axis=1)
profs = {}
for c in range(3):
    m = cl == c
    profs[c] = float(ma[m].mean())
ordr = sorted(profs, key=lambda c: profs[c])
rmap = {ordr[0]: "expansion_competent", ordr[1]: "committed", ordr[2]: "terminal"}
y = np.array([rmap[c] for c in cl])

# 30-gene panel
qc_data = json.loads((PROJ / "p3_qc_panel/output/qc_panel_239.json").read_text())
gp = qc_data.get("panel_genes", qc_data.get("genes", []))
if isinstance(gp[0], dict):
    gp = [g["gene"] if isinstance(g, dict) else g for g in gp]

panel_genes_in_data = [g for g in gp if g in tpm_df.index]
X_panel = np.zeros((len(ch), len(panel_genes_in_data)), dtype=np.float32)
for i, g in enumerate(panel_genes_in_data):
    gi = list(tpm_df.index).index(g)
    X_panel[:, i] = Xr[:, gi]

print(f"Data: {len(ch)} samples, {len(panel_genes_in_data)} panel genes")
print(f"States: {dict(Counter(y))}")
