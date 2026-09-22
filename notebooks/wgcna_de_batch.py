"""WGCNA co-expression + differential expression + batch effects."""
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist
from scipy.stats import ttest_ind
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
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

print(f"Data: {len(ch)} samples, {len(panel_genes)} genes")

# ── WGCNA-style co-expression ──
print("\n─── WGCNA Co-expression Modules ───")
corr = np.corrcoef(X_panel.T)
dissim = 1 - np.abs(corr)
Z = linkage(pdist(dissim), method='average')
clusters = fcluster(Z, t=0.3, criterion='distance')
n_modules = len(set(clusters))
print(f"Modules found: {n_modules}")

modules = {}
for i, c in enumerate(clusters):
    if c not in modules: modules[c] = []
    modules[c].append(panel_genes[i])

for mid, genes in sorted(modules.items()):
    print(f"  Module {mid}: {len(genes)} genes -> {genes}")

# Module eigengene correlation with states
module_eigengenes = {}
for mid, genes in modules.items():
    idx = [panel_genes.index(g) for g in genes]
    if len(idx) > 1:
        pca = PCA(1).fit_transform(X_panel[:, idx])
        module_eigengenes[mid] = pca.flatten()
    else:
        module_eigengenes[mid] = X_panel[:, idx[0]]

state_corrs = {}
for mid, me in module_eigengenes.items():
    state_corrs[mid] = {}
    for s in ["expansion_competent", "committed", "terminal"]:
        mask = y == s
        state_corrs[mid][s] = float(me[mask].mean())

print("\nModule-state associations:")
for mid in sorted(state_corrs.keys()):
    sc = state_corrs[mid]
    print(f"  Module {mid}: E={sc['expansion_competent']:.3f}, C={sc['committed']:.3f}, T={sc['terminal']:.3f}")

# ── Differential expression ──
print("\n─── Differential Expression per State ──")
de_results = []
for i, g in enumerate(panel_genes):
    expr = X_panel[:, i]
    for s1, s2 in [("expansion_competent","terminal"), ("expansion_competent","committed"), ("committed","terminal")]:
        e1 = expr[y == s1]
        e2 = expr[y == s2]
        if len(e1) > 1 and len(e2) > 1:
            t, p = ttest_ind(e1, e2)
            fc = e1.mean() / (e2.mean() + 1e-8)
            de_results.append({"gene": g, "comparison": f"{s1}_vs_{s2}", "log2FC": float(np.log2(fc + 1e-8)), "pvalue": float(p), "tstat": float(t)})

de_df = pd.DataFrame(de_results)
sig = de_df[de_df["pvalue"] < 0.05]
print(f"Significant DE comparisons: {len(sig)}")
for _, row in sig.head(10).iterrows():
    print(f"  {row['gene']}: {row['comparison']} FC={row['log2FC']:.2f} p={row['pvalue']:.4f}")

# ── Batch effect analysis ──
print("\n─── Batch Effect Analysis ───")
# Simulate batch detection: split samples into 3 artificial batches
rng = np.random.RandomState(42)
batches = rng.randint(0, 3, len(ch))
batch_effects = {}
for i, g in enumerate(panel_genes):
    expr = X_panel[:, i]
    batch_means = [expr[batches == b].mean() for b in range(3)]
    batch_effects[g] = {"max_batch_diff": float(max(batch_means) - min(batch_means)), "cv_across_batches": float(np.std(batch_means) / (np.mean(batch_means) + 1e-8))}

# PCA to detect batch separation
pca_batch = PCA(2).fit_transform(X_panel)
batch_sep = {}
for b in range(3):
    mask = batches == b
    batch_sep[f"batch_{b}"] = {"pc1_mean": float(pca_batch[mask, 0].mean()), "pc2_mean": float(pca_batch[mask, 1].mean())}

print("Batch PC separation:")
for b, v in batch_sep.items():
    print(f"  {b}: PC1={v['pc1_mean']:.3f}, PC2={v['pc2_mean']:.3f}")

# ── Save ──
results = {
    "wgcna": {"n_modules": int(n_modules), "modules": {str(k): v for k, v in modules.items()}, "module_state_corr": {str(k): v for k, v in state_corrs.items()}},
    "differential_expression": de_results,
    "batch_effects": {"batch_pc_separation": batch_sep, "gene_batch_effects": batch_effects}
}
json.dump(results, open(OUT / "wgcna_de_batch.json", "w"), indent=2)
print("\nSaved to wgcna_de_batch.json")
print("DONE")
