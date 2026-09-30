"""P2: Multi-Omic State Map for Stem Cell Manufacturing Readiness."""
import json, warnings
from pathlib import Path
import numpy as np, pandas as pd, yaml
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import LeaveOneGroupOut, cross_val_score
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
PROJ = Path(__file__).resolve().parents[1]
CFG = yaml.safe_load((PROJ / "shared_data/data_paths.yaml").read_text())
OUT = PROJ / "p2_state_map/output"
OUT.mkdir(exist_ok=True)

print("=" * 60)
print("P2: MULTI-OMIC STATE MAP")
print("=" * 60)

# 1. Load and align
print("\n1. Loading and aligning RNA + flux...")
norm = pd.read_csv(CFG["rna_seq"]["normalized_counts"], index_col=0)
flux_raw = pd.read_csv(CFG["flux"]["metaflux"], index_col=0)
meta = pd.read_csv(CFG["metadata"]["geo"])

# Build SRR -> sample_name map
srr_to_sample = dict(zip(meta["SRR_Accession"], meta["Sample_ID"]))
sample_to_condition = dict(zip(meta["Sample_ID"], meta["Condition"]))

# Filter flux to only SRRs present in metadata
flux_srr_in_meta = [c for c in flux_raw.columns if c in srr_to_sample]
flux_aligned = flux_raw[flux_srr_in_meta].copy()
flux_aligned.columns = [srr_to_sample[c] for c in flux_aligned.columns]

# Align to RNA sample order
rna_samples = norm.columns.tolist()
common = [s for s in rna_samples if s in flux_aligned.columns]
print(f"   RNA samples: {len(rna_samples)}, Flux samples: {len(flux_aligned.columns)}, Common: {len(common)}")

X_rna = np.log1p(norm[common].values.T)
X_flux = flux_aligned[common].values.T
conditions = np.array([sample_to_condition[s] for s in common])
print(f"   RNA: {X_rna.shape}, Flux: {X_flux.shape}")

# 2. Joint embedding
print("\n2. Joint RNA+flux embedding...")
X_rna_s = StandardScaler().fit_transform(X_rna)
X_flux_s = StandardScaler().fit_transform(X_flux)
pca_rna = PCA(n_components=5).fit_transform(X_rna_s)
pca_flux = PCA(n_components=5).fit_transform(X_flux_s)
X_joint = np.hstack([pca_rna, pca_flux])
print(f"   Joint: {X_joint.shape}")

# 3. Define readiness labels
print("\n3. Defining readiness labels...")
label_map = {}
for c in np.unique(conditions):
    if "TSCM" in c: label_map[c] = "expansion_competent"
    elif "TA" in c: label_map[c] = "committed"
    elif "TUA" in c: label_map[c] = "terminal"
readiness = np.array([label_map[c] for c in conditions])
print(f"   Labels: {np.unique(readiness, return_counts=True)}")

# 4. Cluster and visualize
print("\n4. Clustering...")
km = KMeans(n_clusters=3, random_state=42, n_init=10)
clusters = km.fit_predict(X_joint)
pca_viz = PCA(n_components=2).fit_transform(X_joint)

fig, axes = plt.subplots(1, 2, figsize=(14, 5))
colors = {"expansion_competent": "#2ecc71", "committed": "#f39c12", "terminal": "#e74c3c"}
for label in np.unique(readiness):
    mask = readiness == label
    axes[0].scatter(pca_viz[mask, 0], pca_viz[mask, 1], c=colors[label], label=label, s=120, edgecolors="k")
axes[0].set_title("True Readiness States"); axes[0].legend()
for c in range(3):
    mask = clusters == c
    axes[1].scatter(pca_viz[mask, 0], pca_viz[mask, 1], label=f"Cluster {c}", s=120, edgecolors="k")
axes[1].set_title("K-Means Clusters (k=3)"); axes[1].legend()
plt.tight_layout(); fig.savefig(OUT / "state_map.png", dpi=150)
print("   Saved state_map.png")

# 5. Classifier
print("\n5. Readiness classifier...")
logo = LeaveOneGroupOut()
groups = np.array([c.split("_")[0] for c in conditions])
lr = LogisticRegression(max_iter=1000)
scores = cross_val_score(lr, X_joint, readiness, cv=logo, groups=groups)
print(f"   Leave-one-sex-out accuracy: {scores.mean():.2f} (+/- {scores.std():.2f})")

# 6. Top discriminative features
print("\n6. Top discriminative features...")
lr_full = LogisticRegression(max_iter=1000).fit(X_joint, readiness)
coef_norms = np.linalg.norm(lr_full.coef_, axis=0)
top_idx = np.argsort(coef_norms)[-10:][::-1]
print(f"   Top dims: {top_idx.tolist()}")

# 7. Save
results = {
    "n_samples": int(len(X_joint)),
    "conditions": np.unique(conditions).tolist(),
    "readiness_labels": {c: label_map[c] for c in label_map},
    "leave_one_sex_out_accuracy": float(scores.mean()),
    "accuracy_std": float(scores.std()),
    "top_dims": top_idx.tolist(),
}
with open(OUT / "state_map_results.json", "w") as f:
    json.dump(results, f, indent=2)

print(f"\nOutputs: {OUT}")
print("=" * 60)
print("P2 COMPLETE")
print("=" * 60)