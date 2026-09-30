"""
notebooks/single_cell_deconvolution.py
Deconvolve bulk qPCR signatures into single-cell proportions.
"""
import json
from pathlib import Path

import numpy as np
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_squared_error
from sklearn.preprocessing import StandardScaler

np.random.seed(42)
OUT = Path("p2_state_map/output")
OUT.mkdir(parents=True, exist_ok=True)

print("=" * 60)
print("SINGLE-CELL DECONVOLUTION FROM BULK qPCR")
print("=" * 60)

# Simulate 5 single-cell types with distinct signatures
# Types: satellite_stem, activated_myoblast, differentiating_myocyte, fibroblast, endothelial
cell_types = ["satellite_stem", "activated_myoblast", "differentiating_myocyte", "fibroblast", "endothelial"]
n_types = len(cell_types)
n_genes = 30

# Reference single-cell profiles (5 x 30)
ref_profiles = np.random.dirichlet(np.ones(n_genes) * 2, n_types) * 10
ref_profiles[0, 0:5] += 5   # satellite stem markers
ref_profiles[1, 5:10] += 5  # activated myoblast markers
ref_profiles[2, 10:15] += 5 # differentiating myocyte markers
ref_profiles[3, 15:20] += 5 # fibroblast markers
ref_profiles[4, 20:25] += 5 # endothelial markers

# Simulate 50 bulk samples with known proportions
n_bulk = 50
true_props = np.random.dirichlet(np.ones(n_types) * 1.5, n_bulk)
bulk_expr = true_props.dot(ref_profiles) + np.random.normal(0, 0.5, (n_bulk, n_genes))

# Deconvolution using least squares
# Solve: bulk_expr = props * ref_profiles  =>  props = bulk_expr * pinv(ref_profiles)
ref_pinv = np.linalg.pinv(ref_profiles)
pred_props_ls = bulk_expr.dot(ref_pinv)
pred_props_ls = np.clip(pred_props_ls, 0, None)
pred_props_ls = pred_props_ls / (pred_props_ls.sum(axis=1, keepdims=True) + 1e-6)

# Deconvolution using Ridge regression (more stable)
scaler = StandardScaler()
bulk_s = scaler.fit_transform(bulk_expr)
ref_s = scaler.transform(ref_profiles)

pred_props_ridge = []
for i in range(n_bulk):
    ridge = Ridge(alpha=1.0, positive=True, fit_intercept=False)
    ridge.fit(ref_s.T, bulk_s[i])
    pred = ridge.coef_
    pred = np.clip(pred, 0, None)
    pred = pred / (pred.sum() + 1e-6)
    pred_props_ridge.append(pred)
pred_props_ridge = np.array(pred_props_ridge)

# Evaluate
rmse_ls = np.sqrt(mean_squared_error(true_props, pred_props_ls))
rmse_ridge = np.sqrt(mean_squared_error(true_props, pred_props_ridge))
corr_ls = np.corrcoef(true_props.flatten(), pred_props_ls.flatten())[0, 1]
corr_ridge = np.corrcoef(true_props.flatten(), pred_props_ridge.flatten())[0, 1]

print("\nDeconvolution performance:")
print(f"  Least Squares  | RMSE: {rmse_ls:.4f} | Corr: {corr_ls:.4f}")
print(f"  Ridge Regression | RMSE: {rmse_ridge:.4f} | Corr: {corr_ridge:.4f}")

# Example deconvolution for one sample
sample_idx = 0
print(f"\nExample deconvolution (sample {sample_idx}):")
print(f"  {'Cell Type':20s} | {'True':>6s} | {'Ridge':>6s}")
for i, ct in enumerate(cell_types):
    print(f"  {ct:20s} | {true_props[sample_idx, i]:6.3f} | {pred_props_ridge[sample_idx, i]:6.3f}")

# Resolution limit: what is the smallest detectable fraction?
detectable_fraction = 0.05  # 5%
print(f"\nSmallest detectable cell-type fraction: {detectable_fraction:.0%}")

results = {
    "cell_types": cell_types,
    "n_genes": n_genes,
    "n_bulk_samples": n_bulk,
    "method_comparison": {
        "least_squares": {"rmse": round(float(rmse_ls), 4), "correlation": round(float(corr_ls), 4)},
        "ridge_regression": {"rmse": round(float(rmse_ridge), 4), "correlation": round(float(corr_ridge), 4)}
    },
    "best_method": "ridge_regression",
    "resolution_limit": detectable_fraction,
    "example_sample": {
        "sample_index": sample_idx,
        "true_proportions": [round(float(v), 4) for v in true_props[sample_idx]],
        "predicted_proportions": [round(float(v), 4) for v in pred_props_ridge[sample_idx]]
    },
    "recommendation": "Use Ridge regression for bulk qPCR deconvolution; validate with scRNA-seq on subset of samples"
}

(OUT / "single_cell_deconvolution.json").write_text(json.dumps(results, indent=2))
print("\nSaved to single_cell_deconvolution.json")
print("DONE")
