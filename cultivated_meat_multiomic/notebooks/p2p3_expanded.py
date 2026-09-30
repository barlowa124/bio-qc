"""P2+P3 Expanded: Multi-Omic State Map + QC Panel with 239 aligned samples."""
import json, warnings
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.linear_model import LogisticRegression
from sklearn.feature_selection import SelectFromModel
from sklearn.model_selection import cross_val_score, StratifiedKFold
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
PROJ = Path(__file__).resolve().parents[1]
OUT = PROJ / "p2_state_map/output"
OUT.mkdir(exist_ok=True)
OUT3 = PROJ / "p3_qc_panel/output"
OUT3.mkdir(exist_ok=True)

DATA = Path(r"C:\Users\asdf\CascadeProjects\rao_lab_ml\nn_pipeline\combined_data\metaflux_scfea_pseudobulk")

print("=" * 60)
print("P2+P3 EXPANDED: 239-SAMPLE ANALYSIS")
print("=" * 60)

# 1. Load aligned data
print("\n1. Loading 239-sample aligned data...")
tpm = pd.read_csv(DATA / "scfea_pseudobulk_tpm.csv", index_col=0)
flux = pd.read_csv(DATA / "metaflux_combined_flux_matrix_all.csv", index_col=0)
info = pd.read_csv(DATA / "sample_info.csv")

common = sorted(set(tpm.columns) & set(flux.columns))
print(f"   TPM: {tpm.shape}, Flux: {flux.shape}, Common: {len(common)}")

X_rna = np.log1p(tpm[common].values.T).astype(np.float32)
X_flux = flux[common].values.T.astype(np.float32)
gene_names = tpm.index.values
print(f"   RNA: {X_rna.shape}, Flux: {X_flux.shape}")

# 2. Filter genes
print("\n2. Filtering low-variance genes...")
gvar = X_rna.var(axis=0)
keep = gvar > np.percentile(gvar, 25)
X_rna_f = X_rna[:, keep]
kept = gene_names[keep]
print(f"   Kept {X_rna_f.shape[1]} genes")

# 3. Joint embedding
print("\n3. Joint PCA embedding...")
X_rna_s = StandardScaler().fit_transform(X_rna_f)
X_flux_s = StandardScaler().fit_transform(X_flux)
pca_r = PCA(n_components=15).fit_transform(X_rna_s)
pca_f = PCA(n_components=15).fit_transform(X_flux_s)
X_joint = np.hstack([pca_r, pca_f])
print(f"   Joint: {X_joint.shape}")

# 4. Cluster
print("\n4. Clustering into 3 manufacturing states...")
km = KMeans(n_clusters=3, random_state=42, n_init=10)
clusters = km.fit_predict(X_joint)
meta_act = X_flux_s.mean(axis=1)

profiles = {}
for c in range(3):
    m = clusters == c
    profiles[c] = {"size": int(m.sum()), "metabolic": float(meta_act[m].mean())}

order = sorted(profiles, key=lambda c: profiles[c]["metabolic"])
rmap = {order[0]: "expansion_competent", order[1]: "committed", order[2]: "terminal"}
readiness = np.array([rmap[c] for c in clusters])

for c in range(3):
    print(f"   Cluster {c} ({rmap[c]}): n={profiles[c]['size']}, metabolic={profiles[c]['metabolic']:.3f}")

# 5. Visualize
print("\n5. Visualizing...")
pv = PCA(n_components=2).fit_transform(X_joint)
fig, ax = plt.subplots(1, 2, figsize=(14, 5))
cols = {"expansion_competent": "#2ecc71", "committed": "#f39c12", "terminal": "#e74c3c"}
for lab in np.unique(readiness):
    m = readiness == lab
    ax[0].scatter(pv[m, 0], pv[m, 1], c=cols[lab], label=f"{lab} ({m.sum()})", s=20, alpha=0.5)
ax[0].set_title(f"Manufacturing Readiness (n={len(X_joint)})"); ax[0].legend(markerscale=3)
for c in range(3):
    m = clusters == c
    ax[1].scatter(pv[m, 0], pv[m, 1], label=f"Cluster {c}", s=20, alpha=0.5)
ax[1].set_title("K-Means (k=3)"); ax[1].legend(markerscale=3)
plt.tight_layout(); fig.savefig(OUT / "state_map_239.png", dpi=150)
print("   Saved state_map_239.png")

# 6. Classifier
print("\n6. Readiness classifier...")
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
lr = LogisticRegression(max_iter=2000)
sc = cross_val_score(lr, X_joint, readiness, cv=cv)
print(f"   5-fold CV: {sc.mean():.3f} +/- {sc.std():.3f}")

# 7. QC Panel
print("\n7. QC biomarker panel (L1)...")
lr_l1 = LogisticRegression(penalty="l1", solver="saga", max_iter=5000, C=0.03)
sel = SelectFromModel(lr_l1, max_features=30, threshold=-np.inf)
sel.fit(X_rna_s, readiness)
sm = sel.get_support()
sg = [str(kept[i]) for i in range(len(kept)) if sm[i]]
print(f"   Selected {len(sg)} genes")

Xs = X_rna_s[:, sm]
lf = LogisticRegression(max_iter=2000)
sq = cross_val_score(lf, Xs, readiness, cv=cv)
print(f"   5-fold CV (30-gene): {sq.mean():.3f} +/- {sq.std():.3f}")

lf.fit(Xs, readiness)
print("   Top 10 genes:")
for i in range(min(10, len(sg))):
    cs = " ".join([f"{lf.coef_[j][i]:.2f}" for j in range(len(lf.classes_))])
    print(f"     {sg[i]:20s} [{cs}]")

# 8. Save
json.dump({"n": int(len(X_joint)), "n_genes": int(X_rna_f.shape[1]), "clusters": {str(k): v for k, v in profiles.items()}, "readiness": {str(k): v for k, v in rmap.items()}, "cv_acc": float(sc.mean()), "cv_std": float(sc.std())}, open(OUT / "state_map_239.json", "w"), indent=2)
json.dump({"n_genes": len(sg), "genes": sg, "cv_acc": float(sq.mean()), "cv_std": float(sq.std())}, open(OUT3 / "qc_panel_239.json", "w"), indent=2)

print(f"\nDone. Outputs: {OUT}, {OUT3}")
print("=" * 60)