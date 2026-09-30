"""Build bovine-specific state map and compare with human."""
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

print("=" * 60)
print("BOVINE-SPECIFIC STATE MAP + CROSS-SPECIES COMPARISON")
print("=" * 60)

# ── 1. Load and prep bovine data ──────────────────────────
print("\n1. Preparing bovine data...")
sf = pd.read_csv(CROSS / "GSE173199/sf_diff_counts.csv", index_col=0)
tc = pd.read_csv(CROSS / "GSE173199/timecourse_counts.csv", index_col=0)
bovine = pd.concat([sf, tc], axis=1)

symbol_map = json.loads((CROSS / "bovine_ensembl_to_symbol.json").read_text())
new_index = [symbol_map.get(g.split(".")[0], g) for g in bovine.index.tolist()]
bovine.index = new_index
bovine = bovine[~bovine.index.duplicated(keep='first')]

# CPM + log
bov_cpm = bovine.values / bovine.values.sum(axis=0) * 1e6
bov_log = np.log1p(bov_cpm.T).astype(np.float32)

# Filter low-var genes
gvar = bov_log.var(axis=0)
keep = gvar > np.percentile(gvar, 25)
bov_f = bov_log[:, keep]
bov_genes = bovine.index[keep].tolist()
print(f"   {bov_f.shape[1]} genes after variance filter")

# ── 2. Build bovine state map ─────────────────────────────
print("\n2. Building bovine state map...")
bov_s = StandardScaler().fit_transform(bov_f)
bov_pca = PCA(n_components=15).fit_transform(bov_s)

# K-means with k=2,3,4 to find best separation
for k in [2, 3, 4]:
    km = KMeans(n_clusters=k, random_state=42, n_init=10)
    cl = km.fit_predict(bov_pca)
    inertia = km.inertia_
    counts = Counter(cl)
    print(f"   k={k}: inertia={inertia:.1f}, distribution={dict(counts)}")

# Use k=3 for consistency with human
km = KMeans(n_clusters=3, random_state=42, n_init=10)
bov_clusters = km.fit_predict(bov_pca)

# Characterize clusters by mean expression
bov_meta = bov_s.mean(axis=1)  # mean normalized expression as metabolic proxy
profiles = {}
for c in range(3):
    m = bov_clusters == c
    profiles[c] = {"size": int(m.sum()), "metabolic": float(bov_meta[m].mean())}
order = sorted(profiles, key=lambda c: profiles[c]["metabolic"])
rmap = {order[0]: "expansion_competent", order[1]: "committed", order[2]: "terminal"}
bov_readiness = np.array([rmap[c] for c in bov_clusters])

bov_counts = Counter(bov_readiness)
print("\n   Bovine state distribution:")
for s in ["expansion_competent", "committed", "terminal"]:
    print(f"     {s}: {bov_counts.get(s, 0)} ({bov_counts.get(s,0)/len(bov_readiness):.1%})")

# ── 3. Per-condition breakdown ────────────────────────────
print("\n3. Per-sample states:")
for sample, state in zip(bovine.columns, bov_readiness):
    print(f"     {sample:20s} -> {state}")

print("\n   Per-condition summary:")
conditions = {}
for sample, state in zip(bovine.columns, bov_readiness):
    parts = sample.split("_")
    cond = "_".join(parts[:-1]) if len(parts) > 1 else sample
    if cond not in conditions:
        conditions[cond] = []
    conditions[cond].append(state)

for cond, states in sorted(conditions.items()):
    c = Counter(states)
    print(f"     {cond}: {dict(c)}")

# ── 4. Cross-species: conserved marker analysis ───────────
print("\n4. Cross-species conserved marker analysis...")

conserved = [
    "PAX7", "MYF5", "MYOD1", "MYOG", "MYH3", "DES", "MYL1", "NEB",
    "POU5F1", "SOX2", "NANOG", "MYC", "KLF4",
    "CD44", "CD90", "CD105",
    "GAPDH", "ACTB", "HPRT1",
]

# Get expression of conserved markers in bovine
bov_markers = {}
for gene in conserved:
    if gene in bovine.index:
        expr = bov_log[:, list(bovine.index).index(gene)]
        bov_markers[gene] = {
            "mean": float(expr.mean()),
            "by_state": {s: float(expr[bov_readiness == s].mean()) for s in ["expansion_competent", "committed", "terminal"]}
        }

# Get human marker expression
human_tpm = pd.read_csv(r"C:\Users\asdf\CascadeProjects\rao_lab_ml\nn_pipeline\combined_data\metaflux_scfea_pseudobulk\scfea_pseudobulk_tpm.csv", index_col=0)
flux_df = pd.read_csv(r"C:\Users\asdf\CascadeProjects\rao_lab_ml\nn_pipeline\combined_data\metaflux_scfea_pseudobulk\metaflux_combined_flux_matrix_all.csv", index_col=0)
tpm_df = pd.read_csv(r"C:\Users\asdf\CascadeProjects\rao_lab_ml\nn_pipeline\combined_data\metaflux_scfea_pseudobulk\scfea_pseudobulk_tpm.csv", index_col=0)

ch = sorted(set(tpm_df.columns) & set(flux_df.columns))
Xr = np.log1p(tpm_df[ch].values.T).astype(np.float32)
Xf = flux_df[ch].values.T.astype(np.float32)
gv2 = Xr.var(axis=0)
kp2 = gv2 > np.percentile(gv2, 25)
Xr_f2 = Xr[:, kp2]
Xr_s2 = StandardScaler().fit_transform(Xr_f2)
Xf_s2 = StandardScaler().fit_transform(Xf)

pr2 = PCA(n_components=15).fit_transform(Xr_s2)
pf2 = PCA(n_components=15).fit_transform(Xf_s2)
Xj2 = np.hstack([pr2, pf2])

km2 = KMeans(n_clusters=3, random_state=42, n_init=10)
cl2 = km2.fit_predict(Xj2)
ma2 = Xf_s2.mean(axis=1)
profs2 = {}
for c in range(3):
    m = cl2 == c
    profs2[c] = float(ma2[m].mean())
ordr2 = sorted(profs2, key=lambda c: profs2[c])
rmap2 = {ordr2[0]: "expansion_competent", ordr2[1]: "committed", ordr2[2]: "terminal"}
hum_readiness = np.array([rmap2[c] for c in cl2])

hum_markers = {}
hum_log = np.log1p(tpm_df[ch].values.T).astype(np.float32)
for gene in conserved:
    if gene in tpm_df.index:
        gi = list(tpm_df.index).index(gene)
        expr = hum_log[:, gi]
        hum_markers[gene] = {
            "mean": float(expr.mean()),
            "by_state": {s: float(expr[hum_readiness == s].mean()) for s in ["expansion_competent", "committed", "terminal"]}
        }

# Compare
print("\n   Marker gene comparison (expansion / committed / terminal):")
print(f"   {'Gene':12s} {'Bovine':30s} {'Human':30s}")
print(f"   {'-'*12} {'-'*30} {'-'*30}")
for gene in conserved:
    if gene in bov_markers and gene in hum_markers:
        b = bov_markers[gene]["by_state"]
        h = hum_markers[gene]["by_state"]
        b_str = f"{b['expansion_competent']:.2f}/{b['committed']:.2f}/{b['terminal']:.2f}"
        h_str = f"{h['expansion_competent']:.2f}/{h['committed']:.2f}/{h['terminal']:.2f}"
        print(f"   {gene:12s} {b_str:30s} {h_str:30s}")

# ── 5. Save ───────────────────────────────────────────────
results = {
    "bovine_state_map": {
        "n_samples": len(bov_readiness),
        "distribution": {s: bov_counts.get(s, 0) for s in ["expansion_competent", "committed", "terminal"]},
        "sample_states": dict(zip(bovine.columns.tolist(), [str(s) for s in bov_readiness])),
        "condition_summary": {cond: {k: int(v) for k, v in Counter(states).items()} for cond, states in conditions.items()},
    },
    "cross_species_markers": {
        "bovine": bov_markers,
        "human": hum_markers,
    },
    "common_genes": len([g for g in conserved if g in bov_markers and g in hum_markers]),
}
json.dump(results, open(OUT / "cross_species_comparison.json", "w"), indent=2)
print(f"\nSaved to {OUT}/cross_species_comparison.json")
print("=" * 60)
print("COMPLETE")
print("=" * 60)