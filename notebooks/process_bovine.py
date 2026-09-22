"""Process bovine data v3: use pre-built mapping, fix alignment."""
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
print("BOVINE PROCESSING V3 (using pre-built mapping)")
print("=" * 60)

# ── 1. Load bovine counts + apply pre-built mapping ───────
print("\n1. Loading and mapping bovine genes...")
sf = pd.read_csv(CROSS / "GSE173199/sf_diff_counts.csv", index_col=0)
tc = pd.read_csv(CROSS / "GSE173199/timecourse_counts.csv", index_col=0)
bovine = pd.concat([sf, tc], axis=1)

# Load pre-built mapping
symbol_map = json.loads((CROSS / "bovine_ensembl_to_symbol.json").read_text())
new_index = [symbol_map.get(g.split(".")[0], g) for g in bovine.index.tolist()]
bovine.index = new_index
bovine = bovine[~bovine.index.duplicated(keep='first')]
print(f"   {bovine.shape[0]} genes x {bovine.shape[1]} samples")

# ── 2. Align with human ───────────────────────────────────
print("\n2. Aligning with human genes...")
human_tpm = pd.read_csv(r"C:\Users\asdf\CascadeProjects\rao_lab_ml\nn_pipeline\combined_data\metaflux_scfea_pseudobulk\scfea_pseudobulk_tpm.csv", index_col=0)
human_genes = set(human_tpm.index.tolist())
bovine_symbols = set(bovine.index.tolist())

# Direct match
common = sorted(human_genes & bovine_symbols)
# Case-insensitive
hlower = {g.upper(): g for g in human_genes}
blower = {g.upper(): g for g in bovine_symbols}
cl = sorted(set(hlower.keys()) & set(blower.keys()))
common_ci = [hlower[g] for g in cl]

# Merge both
common_all = sorted(set(common) | set(common_ci))
print(f"   Common genes: {len(common_all)}")

bov_sub = bovine.loc[[g for g in common_all if g in bovine.index]]
hum_sub = human_tpm.loc[[g for g in common_all if g in human_tpm.index]]
print(f"   Bovine: {bov_sub.shape}, Human: {hum_sub.shape}")

# ── 3. Normalize ──────────────────────────────────────────
print("\n3. Normalizing...")
bov_cpm = bov_sub.values / bov_sub.values.sum(axis=0) * 1e6
bov_log = np.log1p(bov_cpm.T).astype(np.float32)
hum_log = np.log1p(hum_sub.values.T).astype(np.float32)

scaler = StandardScaler().fit(hum_log)
hum_s = scaler.transform(hum_log)
bov_s = scaler.transform(bov_log)

pca = PCA(n_components=15).fit(hum_s)
hum_pca = pca.transform(hum_s)
bov_pca = pca.transform(bov_s)
print(f"   Human PCA: {hum_pca.shape}, Bovine PCA: {bov_pca.shape}")

# ── 4. Get human state labels (aligned to TPM samples) ────
print("\n4. Computing human state labels...")
flux_df = pd.read_csv(r"C:\Users\asdf\CascadeProjects\rao_lab_ml\nn_pipeline\combined_data\metaflux_scfea_pseudobulk\metaflux_combined_flux_matrix_all.csv", index_col=0)
tpm_df = pd.read_csv(r"C:\Users\asdf\CascadeProjects\rao_lab_ml\nn_pipeline\combined_data\metaflux_scfea_pseudobulk\scfea_pseudobulk_tpm.csv", index_col=0)

# Align flux and TPM
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
readiness_flux = np.array([rmap[c] for c in cl])

# Map readiness to TPM sample order
# The TPM has 425 samples, flux has 239. Need to map.
tpm_samples = human_tpm.columns.tolist()
flux_samples = ch  # 239 samples

# For TPM-only classification, use nearest-centroid on the 239 aligned samples
# and then apply to bovine
readiness_for_pca = np.full(len(tpm_samples), "unknown", dtype=object)
for i, s in enumerate(tpm_samples):
    if s in flux_samples:
        j = flux_samples.index(s)
        readiness_for_pca[i] = readiness_flux[j]

# Use only labeled samples for centroids
labeled_idx = [i for i in range(len(tpm_samples)) if readiness_for_pca[i] != "unknown"]
print(f"   Labeled human samples: {len(labeled_idx)}")

# ── 5. Classify bovine ────────────────────────────────────
print("\n5. Classifying bovine samples...")
centroids = {}
for s in ["expansion_competent", "committed", "terminal"]:
    idx = [i for i in labeled_idx if readiness_for_pca[i] == s]
    if idx:
        centroids[s] = hum_pca[idx].mean(axis=0)
    else:
        centroids[s] = np.zeros(15)

bov_states = []
bov_dists = []
for i in range(len(bov_pca)):
    dists = {s: float(np.linalg.norm(bov_pca[i] - centroids[s])) for s in centroids}
    best = min(dists, key=dists.get)
    bov_states.append(best)
    bov_dists.append(dists)

counts = Counter(bov_states)
print(f"\n   Bovine state distribution (n={len(bov_states)}):")
for s in ["expansion_competent", "committed", "terminal"]:
    c = counts.get(s, 0)
    print(f"     {s}: {c} ({c/len(bov_states):.1%})")

# Per-sample
print("\n   Per-sample classification:")
for sample, state in zip(bov_sub.columns, bov_states):
    d = bov_dists[len(bov_states) - 1 - list(reversed(bov_sub.columns.tolist())).index(sample)]
    # simpler:
    pass

for i, (sample, state) in enumerate(zip(bov_sub.columns, bov_states)):
    print(f"     {sample:20s} -> {state}")

# Per-condition summary
print("\n   Per-condition summary:")
conditions = {}
for sample, state in zip(bov_sub.columns, bov_states):
    # Extract condition from sample name
    parts = sample.split("_")
    cond = "_".join(parts[:-1]) if len(parts) > 1 else sample
    if cond not in conditions:
        conditions[cond] = []
    conditions[cond].append(state)

for cond, states in sorted(conditions.items()):
    c = Counter(states)
    print(f"     {cond}: {dict(c)}")

# ── 6. Save ───────────────────────────────────────────────
results = {
    "n_bovine_samples": len(bov_states),
    "n_common_genes": len(common_all),
    "state_distribution": {s: counts.get(s, 0) for s in ["expansion_competent", "committed", "terminal"]},
    "sample_states": dict(zip(bov_sub.columns.tolist(), bov_states)),
    "condition_summary": {cond: dict(Counter(states)) for cond, states in conditions.items()},
    "method": "nearest-centroid in human RNA PCA space",
}
json.dump(results, open(OUT / "bovine_cross_species_results.json", "w"), indent=2)
print(f"\nSaved to {OUT}/bovine_cross_species_results.json")
print("=" * 60)
print("COMPLETE")
print("=" * 60)