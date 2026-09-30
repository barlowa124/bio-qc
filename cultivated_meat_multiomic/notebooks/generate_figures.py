"""Generate all 5 publication figures for the manuscript."""
import warnings
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.linear_model import LogisticRegression
from sklearn.feature_selection import SelectFromModel
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import confusion_matrix
import matplotlib; matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import seaborn as sns

warnings.filterwarnings("ignore")
sns.set_style("whitegrid")
plt.rcParams.update({"font.size": 11, "axes.titlesize": 13, "axes.labelsize": 11})

PROJ = Path(__file__).resolve().parents[1]
OUT = PROJ / "p2_state_map/output"
OUT.mkdir(exist_ok=True)
DATA = Path(r"C:\Users\asdf\CascadeProjects\rao_lab_ml\nn_pipeline\combined_data\metaflux_scfea_pseudobulk")

print("=" * 60)
print("GENERATING 5 PUBLICATION FIGURES")
print("=" * 60)

# ── Load + prep ────────────────────────────────────────────
print("\nLoading data...")
tpm = pd.read_csv(DATA / "scfea_pseudobulk_tpm.csv", index_col=0)
flux = pd.read_csv(DATA / "metaflux_combined_flux_matrix_all.csv", index_col=0)
common = sorted(set(tpm.columns) & set(flux.columns))

X_rna = np.log1p(tpm[common].values.T).astype(np.float32)
X_flux = flux[common].values.T.astype(np.float32)
genes = tpm.index.values

gvar = X_rna.var(axis=0)
keep = gvar > np.percentile(gvar, 25)
X_rna_f = X_rna[:, keep]
kept = genes[keep]

X_rna_s = StandardScaler().fit_transform(X_rna_f)
X_flux_s = StandardScaler().fit_transform(X_flux)
pca_r = PCA(n_components=15).fit_transform(X_rna_s)
pca_f = PCA(n_components=15).fit_transform(X_flux_s)
X_joint = np.hstack([pca_r, pca_f])

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

pv = PCA(n_components=2).fit_transform(X_joint)
cols = {"expansion_competent": "#2ecc71", "committed": "#f39c12", "terminal": "#e74c3c"}

# L1 panel
lr_l1 = LogisticRegression(penalty="l1", solver="saga", max_iter=5000, C=0.03)
sel = SelectFromModel(lr_l1, max_features=30, threshold=-np.inf)
sel.fit(X_rna_s, readiness)
sm = sel.get_support()
sg = [str(kept[i]) for i in range(len(kept)) if sm[i]]
Xs = X_rna_s[:, sm]
lf = LogisticRegression(max_iter=2000).fit(Xs, readiness)

# ── Figure 1: State Map ────────────────────────────────────
print("\nFigure 1: State Map...")
fig, ax = plt.subplots(figsize=(8, 7))
for lab in ["expansion_competent", "committed", "terminal"]:
    m = readiness == lab
    ax.scatter(pv[m, 0], pv[m, 1], c=cols[lab], label=lab.replace("_", " ").title(), s=40, alpha=0.6, edgecolors="k", linewidth=0.3)
ax.set_xlabel("PC1"); ax.set_ylabel("PC2")
ax.set_title("Multi-Omic State Map of Stem Cell Manufacturing Readiness\n(239 samples, joint RNA + metabolic flux embedding)")
ax.legend(markerscale=2, frameon=True)
fig.tight_layout(); fig.savefig(OUT / "fig1_state_map.png", dpi=200)
plt.close()

# ── Figure 2: Cluster Characterization ─────────────────────
print("Figure 2: Cluster Characterization...")
fig, axes = plt.subplots(1, 3, figsize=(14, 5))

# 2a: Metabolic activity per state
states = ["expansion_competent", "committed", "terminal"]
metab = [np.mean(meta_act[readiness == s]) for s in states]
axes[0].bar(states, metab, color=[cols[s] for s in states], edgecolor="k")
axes[0].set_title("Mean Metabolic Activity"); axes[0].set_ylabel("Z-score"); axes[0].tick_params(axis="x", rotation=15)

# 2b: Cluster sizes
sizes = [(readiness == s).sum() for s in states]
axes[1].bar(states, sizes, color=[cols[s] for s in states], edgecolor="k")
axes[1].set_title("Samples per State"); axes[1].set_ylabel("Count"); axes[1].tick_params(axis="x", rotation=15)

# 2c: Top gene expression per state
top5 = sg[:5]
for i, s in enumerate(states):
    m = readiness == s
    expr = X_rna_s[m][:, [list(kept).index(g) for g in top5 if g in kept]].mean(axis=0)
    axes[2].bar(np.arange(len(top5)) + i*0.2, expr, width=0.2, color=cols[s], label=s.replace("_"," ").title(), edgecolor="k")
axes[2].set_xticks(np.arange(len(top5)) + 0.2); axes[2].set_xticklabels(top5, rotation=30, ha="right", fontsize=8)
axes[2].set_title("Top 5 Panel Genes per State"); axes[2].legend(fontsize=7)

fig.tight_layout(); fig.savefig(OUT / "fig2_cluster_characterization.png", dpi=200)
plt.close()

# ── Figure 3: QC Panel Performance ─────────────────────────
print("Figure 3: QC Panel Performance...")
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# 3a: Confusion matrix
y_pred = cross_val_predict(lf, Xs, readiness, cv=StratifiedKFold(5, shuffle=True, random_state=42))
cm = confusion_matrix(readiness, y_pred, labels=states)
sns.heatmap(cm, annot=True, fmt="d", xticklabels=[s.replace("_"," ").title() for s in states], yticklabels=[s.replace("_"," ").title() for s in states], cmap="Blues", ax=axes[0])
axes[0].set_title("Confusion Matrix (5-fold CV)"); axes[0].set_ylabel("True"); axes[0].set_xlabel("Predicted")

# 3b: Per-class metrics
from sklearn.metrics import classification_report

cr = classification_report(readiness, y_pred, output_dict=True)
metrics_data = []
for s in states:
    d = cr[s.replace("_", " ").title()] if s.replace("_", " ").title() in cr else cr.get(s, {})
    metrics_data.append({"State": s.replace("_"," ").title(), "Precision": d.get("precision", 0), "Recall": d.get("recall", 0), "F1": d.get("f1-score", 0)})
mdf = pd.DataFrame(metrics_data).set_index("State")
mdf.plot(kind="bar", ax=axes[1], color=["#3498db", "#2ecc71", "#e74c3c"], edgecolor="k")
axes[1].set_title("Per-Class Metrics"); axes[1].set_ylim(0, 1.05); axes[1].legend(loc="lower right")

fig.tight_layout(); fig.savefig(OUT / "fig3_qc_performance.png", dpi=200)
plt.close()

# ── Figure 4: Gene Coefficient Heatmap ─────────────────────
print("Figure 4: Gene Coefficient Heatmap...")
fig, ax = plt.subplots(figsize=(10, 8))
coef_df = pd.DataFrame(lf.coef_, columns=sg, index=[s.replace("_"," ").title() for s in lf.classes_])
# Take top 20 by absolute sum
top20 = coef_df.abs().sum().nlargest(20).index
sns.heatmap(coef_df[top20], annot=True, fmt=".2f", cmap="RdBu_r", center=0, ax=ax, cbar_kws={"label": "Coefficient"})
ax.set_title("L1-Selected Gene Coefficients per Manufacturing State")
fig.tight_layout(); fig.savefig(OUT / "fig4_gene_heatmap.png", dpi=200)
plt.close()

# ── Figure 5: Study Design Schematic ───────────────────────
print("Figure 5: Study Design Schematic...")
fig, ax = plt.subplots(figsize=(12, 4))
ax.set_xlim(0, 12); ax.set_ylim(0, 4); ax.axis("off")

boxes = [
    (0.5, 2, 2, 1.5, "RNA-seq\n(23,682 genes)", "#3498db"),
    (3, 2, 2, 1.5, "METAFlux\n(13,082 reactions)", "#e74c3c"),
    (5.5, 2, 2, 1.5, "Joint PCA\nEmbedding\n(30 dims)", "#9b59b6"),
    (8, 2.5, 2, 1, "K-Means\n(k=3)", "#f39c12"),
    (8, 1, 2, 1, "L1 Logistic\nRegression", "#2ecc71"),
    (10.5, 2.5, 1.2, 1, "State\nMap", "#f39c12"),
    (10.5, 1, 1.2, 1, "30-Gene\nQC Panel", "#2ecc71"),
]

for x, y, w, h, label, color in boxes:
    rect = mpatches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.1", facecolor=color, edgecolor="k", alpha=0.85)
    ax.add_patch(rect)
    ax.text(x + w/2, y + h/2, label, ha="center", va="center", fontsize=9, color="white", fontweight="bold")

# Arrows
for x1, x2, y in [(2.5, 3, 2.75), (5, 5.5, 2.75), (7.5, 8, 3), (7.5, 8, 1.5), (10, 10.5, 3), (10, 10.5, 1.5)]:
    ax.annotate("", xy=(x2, y), xytext=(x1, y), arrowprops=dict(arrowstyle="->", lw=2, color="#555"))

ax.set_title("Study Design: From Multi-Omic Data to Manufacturing QC Panel", fontsize=14, fontweight="bold", pad=20)
fig.tight_layout(); fig.savefig(OUT / "fig5_study_design.png", dpi=200)
plt.close()

print(f"\nAll 5 figures saved to: {OUT}")
print("=" * 60)
print("FIGURES COMPLETE")
print("=" * 60)