"""Rebuild Y from METAFlux pathway activities + retrain per-domain Ridge."""
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from sklearn.linear_model import RidgeCV
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")
PROJ = Path(__file__).resolve().parents[1]
OUT = PROJ / "p2_state_map/output"
OUT.mkdir(exist_ok=True)

print("=" * 60)
print("METAFLUX Y REBUILD + RIDGE RETRAIN")
print("=" * 60)

# Load pathway-aggregated flux (already mapped!)
pathway_flux = pd.read_csv(r"C:\Users\asdf\CascadeProjects\rao_lab_ml\nn_pipeline\combined_data\metaflux_reaggregated_cleaned.csv", index_col=0)
print(f"\n1. Pathway flux: {pathway_flux.shape}")
print(f"   Pathways: {pathway_flux.columns.tolist()}")

# Load TPM
tpm = pd.read_csv(r"C:\Users\asdf\CascadeProjects\rao_lab_ml\nn_pipeline\combined_data\metaflux_scfea_pseudobulk\scfea_pseudobulk_tpm.csv", index_col=0)
info = pd.read_csv(r"C:\Users\asdf\CascadeProjects\rao_lab_ml\nn_pipeline\combined_data\metaflux_scfea_pseudobulk\sample_info.csv")

# Align
common = sorted(set(tpm.columns) & set(pathway_flux.index))
print(f"   Common samples: {len(common)}")

X = np.log1p(tpm[common].values.T).astype(np.float32)
Y_raw = pathway_flux.loc[common].values.astype(np.float32)

# Remove "Other" column
pathway_cols = pathway_flux.columns.tolist()
keep_cols = [c for c in pathway_cols if c != "Other"]
Y = pathway_flux.loc[common, keep_cols].values.astype(np.float32)
print(f"   X: {X.shape}, Y: {Y.shape} (10 pathways)")

# Filter low-var genes
gvar = X.var(axis=0)
keep = gvar > np.percentile(gvar, 25)
X_f = X[:, keep]
print(f"   After gene filter: {X_f.shape[1]} genes")

# Per-dataset split
datasets = [info[info["sample_id"] == s]["dataset"].values[0] if len(info[info["sample_id"] == s]) > 0 else "unknown" for s in common]
uniq_ds = sorted(set(datasets))

print("\n2. Per-domain Ridge results:")
results = {}
all_y_true = []
all_y_pred = []

for ds in uniq_ds:
    idx = [i for i, d in enumerate(datasets) if d == ds]
    if len(idx) < 5:
        continue
    
    X_ds = X_f[idx]
    Y_ds = Y[idx]
    
    # Standardize
    X_s = StandardScaler().fit_transform(X_ds)
    
    # Per-pathway Ridge
    pathway_corrs = {}
    for j, pw in enumerate(keep_cols):
        y = Y_ds[:, j]
        if y.std() < 1e-8:
            pathway_corrs[pw] = 0.0
            continue
        model = RidgeCV(alphas=[0.1, 1.0, 10.0, 100.0])
        model.fit(X_s, y)
        y_pred = model.predict(X_s)
        corr, _ = pearsonr(y, y_pred)
        pathway_corrs[pw] = float(corr)
    
    mean_corr = np.mean(list(pathway_corrs.values()))
    results[ds] = {"n": len(idx), "mean_corr": mean_corr, "pathways": pathway_corrs}
    print(f"   {ds}: n={len(idx)}, mean_corr={mean_corr:.3f}")
    for pw, c in sorted(pathway_corrs.items(), key=lambda x: -x[1])[:5]:
        print(f"     {pw}: {c:.3f}")

# Overall
overall_mean = np.mean([r["mean_corr"] for r in results.values()])
print(f"\n   Overall mean pathway correlation: {overall_mean:.3f}")

json.dump({"per_domain": results, "overall_mean_corr": overall_mean, "n_samples": len(common), "n_pathways": len(keep_cols), "pathways": keep_cols}, open(OUT / "metaflux_ridge_results.json", "w"), indent=2)
print(f"\nSaved to {OUT}/metaflux_ridge_results.json")
print("=" * 60)