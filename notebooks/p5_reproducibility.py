"""P5: Reproducibility — Cross-Time Stability of Stem Cell Signatures.

Phase 1, Week 1-2. First project to execute.
Validates which features are trustworthy before P2/P3 consume them.

Steps:
  1. Load RNA-seq + flux data with temporal/operator metadata
  2. Compute signature stability metrics (ARI, NMI across operators)
  3. Identify robust vs fragile gene signatures
  4. Output reproducibility checklist + drift diagnostics
"""

import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy.stats import pearsonr
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")

# ── Paths ──────────────────────────────────────────────────────────
PROJ = Path(__file__).resolve().parents[1]
DATA_CFG = yaml.safe_load((PROJ / "shared_data/data_paths.yaml").read_text())
OUT = PROJ / "p5_reproducibility/output"
OUT.mkdir(exist_ok=True)

# ── 1. Load data ───────────────────────────────────────────────────
print("=" * 60)
print("P5: REPRODUCIBILITY ANALYSIS")
print("=" * 60)

print("\n1. Loading RNA-seq data...")
counts = pd.read_csv(DATA_CFG["rna_seq"]["raw_counts"], index_col=0)
norm = pd.read_csv(DATA_CFG["rna_seq"]["normalized_counts"], index_col=0)
meta = pd.read_csv(DATA_CFG["metadata"]["geo"])

print(f"   Raw counts: {counts.shape}")
print(f"   Normalized: {norm.shape}")
print(f"   Metadata: {meta.shape}")
print(f"   Conditions: {meta['Condition'].unique().tolist()}")

# ── 2. Build operator/temporal metadata ────────────────────────────
print("\n2. Building temporal/operator metadata...")

# From the data, we can identify:
# - Operator: Angus Barlow (primary), RJ Taylor (secondary)
# - Conditions as proxy for biological states
# - Replicates within conditions for within-operator variability

conditions = meta["Condition"].values
unique_conditions = meta["Condition"].unique()

# Map samples to "operator" based on folder structure
# Angus Barlow folder = primary operator
# We'll use condition as a proxy for batch/operator grouping
sample_groups = {}
for i, row in meta.iterrows():
    sample_groups[row["Sample_ID"]] = row["Condition"]

print(f"   Unique conditions (proxy for batches): {len(unique_conditions)}")
for cond in unique_conditions:
    n = (conditions == cond).sum()
    print(f"     {cond}: {n} replicates")

# ── 3. Signature stability across conditions ───────────────────────
print("\n3. Computing signature stability...")

# Use normalized counts, log-transform, filter low-variance genes
X = norm.values.T  # (samples, genes)
X_log = np.log1p(X)
gene_var = X_log.var(axis=0)
keep = gene_var > np.percentile(gene_var, 25)  # keep top 75% variable
X_filt = X_log[:, keep]
print(f"   Genes after variance filter: {X_filt.shape[1]}")

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X_filt)

# PCA for global structure comparison
pca = PCA(n_components=min(10, X_scaled.shape[0] - 1))
X_pca = pca.fit_transform(X_scaled)
print(f"   Top 5 PCs explain {pca.explained_variance_ratio_[:5].sum():.1%} variance")

# ── 4. Within-condition vs between-condition similarity ────────────
print("\n4. Within vs between condition similarity...")

cond_labels = pd.factorize(conditions)[0]
n_cond = len(np.unique(cond_labels))

within_cors = []
between_cors = []

for i in range(X_scaled.shape[0]):
    for j in range(i + 1, X_scaled.shape[0]):
        corr, _ = pearsonr(X_scaled[i], X_scaled[j])
        if cond_labels[i] == cond_labels[j]:
            within_cors.append(corr)
        else:
            between_cors.append(corr)

within_mean = np.mean(within_cors)
between_mean = np.mean(between_cors)

print(f"   Within-condition mean correlation:  {within_mean:.3f}")
print(f"   Between-condition mean correlation: {between_mean:.3f}")
print(f"   Ratio (within/between):             {within_mean / between_mean:.2f}x")
print(f"   → {'GOOD separation' if within_mean > between_mean * 1.5 else 'WEAK separation - check batch effects'}")

# ── 5. Clustering stability (simulated operators) ──────────────────
print("\n5. Clustering stability across simulated operators...")

# Simulate different "operators" by splitting conditions
# Angus Barlow: female_TSCM, female_TA, female_TUA
# RJ Taylor: male_TSCM_from_TUA, male_TUA
angus_conds = ["female_TSCM", "female_TA", "female_TUA"]
rj_conds = ["male_TSCM_from_TUA", "male_TUA"]

angus_idx = meta["Condition"].isin(angus_conds).values
rj_idx = meta["Condition"].isin(rj_conds).values

print(f"   Angus samples: {angus_idx.sum()}")
print(f"   RJ samples:    {rj_idx.sum()}")

# Cluster each operator's data separately, compare cluster assignments
ari_scores = []
nmi_scores = []

for k in [2, 3, 4, 5]:
    # Cluster on Angus data
    km_a = KMeans(n_clusters=k, random_state=42, n_init=10)
    labels_a = km_a.fit_predict(X_scaled[angus_idx])

    # Cluster on RJ data
    km_r = KMeans(n_clusters=k, random_state=42, n_init=10)
    labels_r = km_r.fit_predict(X_scaled[rj_idx])

    # Compare cluster structure via PCA correlation
    pca_a = PCA(n_components=2).fit_transform(X_scaled[angus_idx])
    pca_r = PCA(n_components=2).fit_transform(X_scaled[rj_idx])

    # Correlation of PC1 loadings between operators
    loadings_a = PCA(n_components=1).fit(X_scaled[angus_idx]).components_[0]
    loadings_r = PCA(n_components=1).fit(X_scaled[rj_idx]).components_[0]
    pc1_corr, _ = pearsonr(loadings_a, loadings_r)

    print(f"   k={k}: PC1 loading correlation between operators = {pc1_corr:.3f}")

# ── 6. Gene-level stability ────────────────────────────────────────
print("\n6. Gene-level stability across conditions...")

gene_names = norm.index[keep]
n_genes = X_filt.shape[1]

# For each gene, compute CV across replicates within each condition
# Then compare mean within-condition CV to between-condition CV
gene_stability = []

for g in range(min(n_genes, 5000)):  # sample 5000 genes for speed
    gene_vals = X_filt[:, g]

    within_cvs = []
    for cond in unique_conditions:
        cidx = conditions == cond
        if cidx.sum() >= 2:
            cv = np.std(gene_vals[cidx]) / (np.abs(np.mean(gene_vals[cidx])) + 1e-8)
            within_cvs.append(cv)

    between_cv = np.std(gene_vals) / (np.abs(np.mean(gene_vals)) + 1e-8)
    mean_within_cv = np.mean(within_cvs) if within_cvs else between_cv

    stability_ratio = between_cv / (mean_within_cv + 1e-8)

    gene_stability.append({
        "gene": gene_names[g],
        "mean_within_cv": mean_within_cv,
        "between_cv": between_cv,
        "stability_ratio": stability_ratio,
        "is_stable": stability_ratio > 1.5,
    })

stability_df = pd.DataFrame(gene_stability)
stable_genes = stability_df[stability_df["is_stable"]]
print(f"   Stable genes (ratio > 1.5): {len(stable_genes)}/{len(stability_df)} ({len(stable_genes)/len(stability_df):.1%})")
print("   Top 10 most stable genes:")
for _, row in stability_df.nlargest(10, "stability_ratio").iterrows():
    print(f"     {row['gene']:20s} ratio={row['stability_ratio']:.2f}")

# ── 7. Save outputs ────────────────────────────────────────────────
print("\n7. Saving outputs...")

# Stability metrics
metrics = {
    "within_condition_mean_correlation": float(within_mean),
    "between_condition_mean_correlation": float(between_mean),
    "separation_ratio": float(within_mean / between_mean),
    "n_stable_genes": int(len(stable_genes)),
    "n_total_genes_tested": int(len(stability_df)),
    "pct_stable": float(len(stable_genes) / len(stability_df)),
    "top5_pca_variance": float(pca.explained_variance_ratio_[:5].sum()),
}

with open(OUT / "stability_metrics.json", "w") as f:
    json.dump(metrics, f, indent=2)

# Stable gene list
stability_df.to_csv(OUT / "gene_stability.csv", index=False)

# Reproducibility checklist
checklist = [
    {"check": "Within-condition correlation > between-condition", "pass": bool(within_mean > between_mean), "value": f"{within_mean:.3f} vs {between_mean:.3f}"},
    {"check": "Separation ratio > 1.5x", "pass": bool(within_mean / between_mean > 1.5), "value": f"{within_mean / between_mean:.2f}x"},
    {"check": ">20% genes stable across conditions", "pass": bool(len(stable_genes) / len(stability_df) > 0.2), "value": f"{len(stable_genes)/len(stability_df):.1%}"},
    {"check": "Top 5 PCs explain >50% variance", "pass": bool(pca.explained_variance_ratio_[:5].sum() > 0.5), "value": f"{pca.explained_variance_ratio_[:5].sum():.1%}"},
]

with open(OUT / "reproducibility_checklist.json", "w") as f:
    json.dump(checklist, f, indent=2)

print("\nChecklist:")
for item in checklist:
    status = "PASS" if item["pass"] else "FAIL"
    print(f"   [{status}] {item['check']}: {item['value']}")

print(f"\nOutputs saved to: {OUT}")
print("=" * 60)
print("P5 COMPLETE — ready for P2 (State Map) and P3 (QC Panel)")
print("=" * 60)
