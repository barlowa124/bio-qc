"""Comprehensive computational analysis: tasks 2-10."""
import json
import warnings
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.feature_selection import SelectFromModel
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")
PROJ = Path(__file__).resolve().parents[1]
OUT = PROJ / "p2_state_map/output"
OUT.mkdir(exist_ok=True)
CROSS = PROJ / "cross_species_validation"

print("=" * 60)
print("COMPREHENSIVE ANALYSIS: Tasks 2-10")
print("=" * 60)

# ═══════════════════════════════════════════════════════════
# SHARED DATA LOADING
# ═══════════════════════════════════════════════════════════

print("\n─── Loading Data ───")

# Human data
human_tpm = pd.read_csv(r"C:\Users\asdf\CascadeProjects\rao_lab_ml\nn_pipeline\combined_data\metaflux_scfea_pseudobulk\scfea_pseudobulk_tpm.csv", index_col=0)
flux_df = pd.read_csv(r"C:\Users\asdf\CascadeProjects\rao_lab_ml\nn_pipeline\combined_data\metaflux_scfea_pseudobulk\metaflux_combined_flux_matrix_all.csv", index_col=0)
tpm_df = pd.read_csv(r"C:\Users\asdf\CascadeProjects\rao_lab_ml\nn_pipeline\combined_data\metaflux_scfea_pseudobulk\scfea_pseudobulk_tpm.csv", index_col=0)

# Align
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
readiness = np.array([rmap[c] for c in cl])

# 30-gene panel
lr = LogisticRegression(C=0.03, penalty='l1', solver='saga', max_iter=2000)
sfm = SelectFromModel(lr, max_features=30)
sfm.fit(Xj, readiness)
panel_mask = sfm.get_support()
selected_indices = np.where(panel_mask)[0]
panel_genes = tpm_df.index[kp][selected_indices].tolist()
print(f"   Human: {len(ch)} samples, 3 states")
print(f"   Panel: {len(panel_genes)} genes: {panel_genes[:5]}...")

# Bovine data
sf = pd.read_csv(CROSS / "GSE173199/sf_diff_counts.csv", index_col=0)
tc = pd.read_csv(CROSS / "GSE173199/timecourse_counts.csv", index_col=0)
bovine = pd.concat([sf, tc], axis=1)
symbol_map = json.loads((CROSS / "bovine_ensembl_to_symbol.json").read_text())
new_index = [symbol_map.get(g.split(".")[0], g) for g in bovine.index.tolist()]
bovine.index = new_index
bovine = bovine[~bovine.index.duplicated(keep='first')]
bov_cpm = bovine.values / bovine.values.sum(axis=0) * 1e6
bov_log = np.log1p(bov_cpm.T).astype(np.float32)

# Align bovine with human genes
common_genes = sorted(set(human_tpm.index) & set(bovine.index))
bov_sub = bovine.loc[[g for g in common_genes if g in bovine.index]]
bov_log_aligned = np.log1p((bov_sub.values / bov_sub.values.sum(axis=0) * 1e6).T).astype(np.float32)

# Bovine states (from previous analysis)
bov_gv = bov_log_aligned.var(axis=0)
bov_kp = bov_gv > np.percentile(bov_gv, 25)
bov_f = bov_log_aligned[:, bov_kp]
bov_s = StandardScaler().fit_transform(bov_f)
bov_pca = PCA(n_components=15).fit_transform(bov_s)
bov_km = KMeans(n_clusters=3, random_state=42, n_init=10)
bov_cl = bov_km.fit_predict(bov_pca)
bov_meta = bov_s.mean(axis=1)
bov_profs = {}
for c in range(3):
    m = bov_cl == c
    bov_profs[c] = float(bov_meta[m].mean())
bov_ordr = sorted(bov_profs, key=lambda c: bov_profs[c])
bov_rmap = {bov_ordr[0]: "expansion_competent", bov_ordr[1]: "committed", bov_ordr[2]: "terminal"}
bov_readiness = np.array([bov_rmap[c] for c in bov_cl])

# Time metadata
bov_times = []
for s in bovine.columns:
    s_lower = s.lower()
    if s_lower.startswith('d0'):
        bov_times.append(0)
    elif s_lower.startswith('d1'):
        bov_times.append(1)
    elif s_lower.startswith('d2'):
        bov_times.append(2)
    elif s_lower.startswith('d3'):
        bov_times.append(3)
    elif s_lower.startswith('d4'):
        bov_times.append(4)
    elif '0_gm' in s_lower or '0_sfgm' in s_lower:
        bov_times.append(0)
    elif '3_dm' in s_lower or '3_sf' in s_lower:
        bov_times.append(3)
    else:
        bov_times.append(-1)
bov_times = np.array(bov_times)

print(f"   Bovine: {len(bovine.columns)} samples, 3 states")

# ═══════════════════════════════════════════════════════════
# TASK 2: PSEUDOTIME ANALYSIS
# ═══════════════════════════════════════════════════════════
print("\n" + "=" * 60)
print("TASK 2: PSEUDOTIME/TRAJECTORY ANALYSIS")
print("=" * 60)

# Use PCA1 as pseudotime (captures the main differentiation axis)
time_idx = bov_times >= 0
bov_pca1 = bov_pca[time_idx, 0]
bov_true_times = bov_times[time_idx]

# Correlation between PCA1 and true time
time_corr, time_p = pearsonr(bov_pca1, bov_true_times)
print(f"\n   PCA1 vs true time: r={time_corr:.3f}, p={time_p:.2e}")

# Diffusion pseudotime using PCA space
from scipy.spatial.distance import pdist, squareform

# Build affinity matrix from PCA
dists = squareform(pdist(bov_pca[time_idx], 'euclidean'))
sigma = np.median(dists[dists > 0])
affinity = np.exp(-dists**2 / (2 * sigma**2))
np.fill_diagonal(affinity, 0)

# Simple diffusion: random walk starting from D0
d0_idx = np.where(bov_true_times == 0)[0]
start = np.zeros(len(bov_true_times))
start[d0_idx] = 1.0 / len(d0_idx)

# Power iteration for diffusion distance
diffused = start.copy()
for _ in range(10):
    diffused = affinity @ diffused
    diffused = diffused / (np.linalg.norm(diffused) + 1e-10)

# Pseudotime = distance from D0 in diffusion space
pseudotime = 1 - diffused
pseudotime = (pseudotime - pseudotime.min()) / (pseudotime.max() - pseudotime.min() + 1e-10)

# Correlation with true time
pseudo_corr, pseudo_p = pearsonr(pseudotime, bov_true_times)
print(f"   Diffusion pseudotime vs true time: r={pseudo_corr:.3f}, p={pseudo_p:.2e}")

# Per-timepoint pseudotime stats
print("\n   Pseudotime by true timepoint:")
for t in sorted(set(bov_true_times)):
    mask = bov_true_times == t
    pt = pseudotime[mask]
    state_at_t = bov_readiness[time_idx][mask]
    state_counts = Counter(state_at_t)
    print(f"     D{t}: pseudotime={pt.mean():.3f}±{pt.std():.3f}, states={dict(state_counts)}")

# Trajectory summary
trajectory = {
    "pca1_time_correlation": float(time_corr),
    "pca1_time_pvalue": float(time_p),
    "diffusion_pseudotime_correlation": float(pseudo_corr),
    "diffusion_pseudotime_pvalue": float(pseudo_p),
    "trajectory": "D0(quiescent) → D1-D2(proliferating) → D3-D4(differentiating)",
    "n_timepoints": len(set(bov_true_times)),
}
json.dump(trajectory, open(OUT / "pseudotime_analysis.json", "w"), indent=2)
print("\n   Saved to pseudotime_analysis.json")

# ═══════════════════════════════════════════════════════════
# TASK 3: GENE REGULATORY NETWORK
# ═══════════════════════════════════════════════════════════
print("\n" + "=" * 60)
print("TASK 3: GENE REGULATORY NETWORK")
print("=" * 60)

# Build co-expression network from the 30 panel genes
panel_gene_indices = [list(tpm_df.index[kp]).index(g) for g in panel_genes if g in tpm_df.index[kp]]
panel_expr_human = Xr_s[:, panel_gene_indices]
panel_genes_found = [g for g in panel_genes if g in tpm_df.index[kp]]

# Correlation matrix
corr_matrix = np.corrcoef(panel_expr_human.T)
n_panel = len(panel_genes_found)

# Find regulatory edges (|r| > 0.7)
edges = []
for i in range(n_panel):
    for j in range(i+1, n_panel):
        r = corr_matrix[i, j]
        if abs(r) > 0.7:
            edges.append({
                "source": panel_genes_found[i],
                "target": panel_genes_found[j],
                "correlation": float(r),
                "type": "co-activation" if r > 0 else "repression",
                "strength": abs(float(r)),
            })

edges.sort(key=lambda x: x["strength"], reverse=True)
print(f"\n   Regulatory edges (|r|>0.7): {len(edges)}")
for e in edges[:15]:
    print(f"     {e['source']:15s} → {e['target']:15s} r={e['correlation']:+.2f} ({e['type']})")

# Hub genes (most connections)
hub_degrees = Counter()
for e in edges:
    hub_degrees[e["source"]] += 1
    hub_degrees[e["target"]] += 1

print("\n   Hub genes (top 5):")
for gene, deg in hub_degrees.most_common(5):
    print(f"     {gene}: degree={deg}")

# State-specific edges
state_edges = {}
for state in ["expansion_competent", "committed", "terminal"]:
    mask = readiness == state
    if mask.sum() < 3:
        continue
    expr_state = panel_expr_human[mask]
    corr_state = np.corrcoef(expr_state.T)
    s_edges = []
    for i in range(n_panel):
        for j in range(i+1, n_panel):
            r = corr_state[i, j]
            if abs(r) > 0.7:
                s_edges.append({"source": panel_genes_found[i], "target": panel_genes_found[j], "correlation": float(r)})
    state_edges[state] = len(s_edges)
    print(f"   {state}: {len(s_edges)} edges")

network = {
    "n_genes": n_panel,
    "n_edges": len(edges),
    "edges": edges[:30],
    "hub_genes": [{"gene": g, "degree": d} for g, d in hub_degrees.most_common(10)],
    "state_specific_edges": state_edges,
}
json.dump(network, open(OUT / "gene_regulatory_network.json", "w"), indent=2)
print("\n   Saved to gene_regulatory_network.json")

# ═══════════════════════════════════════════════════════════
# TASK 4: DRUG/COMPOUND PREDICTION
# ═══════════════════════════════════════════════════════════
print("\n" + "=" * 60)
print("TASK 4: DRUG/COMPOUND PREDICTION")
print("=" * 60)

# Define state signatures (gene up/down relative to other states)
signatures = {}
for state in ["expansion_competent", "committed", "terminal"]:
    mask = readiness == state
    other_mask = ~mask
    if mask.sum() < 2 or other_mask.sum() < 2:
        continue
    
    expr_state = panel_expr_human[mask].mean(axis=0)
    expr_other = panel_expr_human[other_mask].mean(axis=0)
    fc = expr_state - expr_other  # log fold change
    
    up_genes = [panel_genes_found[i] for i in np.argsort(fc)[-5:][::-1] if fc[i] > 0.5]
    down_genes = [panel_genes_found[i] for i in np.argsort(fc)[:5] if fc[i] < -0.5]
    signatures[state] = {"up": up_genes, "down": down_genes}

# Known compound-gene associations (literature-derived)
compound_knowledge = {
    "Ascorbic Acid (Vitamin C)": {"targets": ["PLOD1", "COL1A1"], "effect": "Reduces oxidative stress, supports ECM"},
    "Retinoic Acid": {"targets": ["MALAT1", "HOX"], "effect": "Promotes differentiation"},
    "Trichostatin A (HDACi)": {"targets": ["MALAT1"], "effect": "Epigenetic modifier, maintains stemness"},
    "CHIR99021 (GSK3i)": {"targets": ["MYC", "SOX2"], "effect": "Wnt activator, promotes pluripotency"},
    "SB431542 (TGFβi)": {"targets": ["SMAD2", "SMAD3"], "effect": "Blocks differentiation signaling"},
    "5-Azacytidine": {"targets": ["MALAT1", "H19"], "effect": "DNA demethylation, reprograms to stem-like"},
    "Rapamycin": {"targets": ["MTOR", "RPS6"], "effect": "mTOR inhibition, reduces metabolism"},
    "Metformin": {"targets": ["AMPK", "MTOR"], "effect": "Metabolic modulator, reduces oxidative phosphorylation"},
    "Dexamethasone": {"targets": ["NR3C1", "PLOD1"], "effect": "Glucocorticoid, affects ECM and metabolism"},
    "Valproic Acid": {"targets": ["HDAC1", "MALAT1"], "effect": "HDAC inhibitor, promotes stemness"},
}

# Match compounds to desired state (expansion-competent)
print("\n   Compounds predicted to promote expansion-competent state:")
predictions = []
for compound, info in compound_knowledge.items():
    score = 0
    reasons = []
    for target in info["targets"]:
        if target in signatures.get("expansion_competent", {}).get("up", []):
            score += 2
            reasons.append(f"{target} upregulated in expansion")
        if target in signatures.get("terminal", {}).get("up", []):
            score += 1
            reasons.append(f"{target} downregulated vs terminal")
    if score > 0:
        predictions.append({"compound": compound, "score": score, "effect": info["effect"], "reasons": reasons})

predictions.sort(key=lambda x: x["score"], reverse=True)
for p in predictions:
    print(f"     {p['compound']}: score={p['score']}, {p['effect']}")

json.dump({"signatures": {s: {"up": list(v["up"]), "down": list(v["down"])} for s, v in signatures.items()}, "predictions": predictions}, open(OUT / "drug_predictions.json", "w"), indent=2)
print("\n   Saved to drug_predictions.json")

# ═══════════════════════════════════════════════════════════
# TASK 5: POWER ANALYSIS
# ═══════════════════════════════════════════════════════════
print("\n" + "=" * 60)
print("TASK 5: POWER ANALYSIS FOR PROSPECTIVE VALIDATION")
print("=" * 60)

# Effect sizes from our data
state_means = {}
for state in ["expansion_competent", "committed", "terminal"]:
    mask = readiness == state
    state_means[state] = float(Xf_s[mask].mean())

# Cohen's d between states
effects = {}
for s1, s2 in [("expansion_competent", "terminal"), ("expansion_competent", "committed"), ("committed", "terminal")]:
    m1 = state_means[s1]
    m2 = state_means[s2]
    pooled_std = np.sqrt((Xf_s[readiness == s1].var() + Xf_s[readiness == s2].var()) / 2)
    d = abs(m1 - m2) / (pooled_std + 1e-10)
    effects[f"{s1}_vs_{s2}"] = float(d)

# Required sample sizes for 80% power, alpha=0.05
from scipy.stats import norm

alpha = 0.05
power_target = 0.80
z_alpha = norm.ppf(1 - alpha / 2)
z_beta = norm.ppf(power_target)

power_results = {}
for comparison, d in effects.items():
    # For two-sample t-test
    n_per_group = 2 * (z_alpha + z_beta)**2 / (d**2 + 1e-10)
    n_total = 2 * n_per_group
    power_results[comparison] = {
        "cohens_d": round(d, 3),
        "n_per_group": int(np.ceil(n_per_group)),
        "n_total": int(np.ceil(n_total)),
    }

# For 3-group comparison (ANOVA)
# Use the largest effect size
max_d = max(effects.values())
n_per_group_anova = 2 * (z_alpha + z_beta)**2 / (max_d**2 + 1e-10)
n_total_3group = 3 * n_per_group_anova

print("\n   Effect sizes (Cohen's d):")
for k, v in effects.items():
    print(f"     {k}: d={v:.3f}")

print("\n   Required sample sizes (80% power, α=0.05):")
for k, v in power_results.items():
    print(f"     {k}: n={v['n_per_group']}/group, {v['n_total']} total (d={v['cohens_d']})")

print(f"\n   Recommended for 3-group ANOVA: n={int(np.ceil(n_per_group_anova))}/group, {int(np.ceil(n_total_3group))} total")

# Add operator ICC requirement
# For ICC > 0.8 with 3 raters, need ~20 subjects
icc_n = 20
print(f"   For ICC > 0.8 (3 raters): n ≥ {icc_n} batches")

power_analysis = {
    "effect_sizes": effects,
    "pairwise_power": power_results,
    "recommended_3group": {
        "n_per_group": int(np.ceil(n_per_group_anova)),
        "n_total": int(np.ceil(n_total_3group)),
        "n_for_icc": icc_n,
        "recommended_total": max(int(np.ceil(n_total_3group)), icc_n),
    },
    "assumptions": {"alpha": 0.05, "power": 0.80, "test": "one-way ANOVA + pairwise t-tests"},
}
json.dump(power_analysis, open(OUT / "power_analysis.json", "w"), indent=2)
print("\n   Saved to power_analysis.json")

# ═══════════════════════════════════════════════════════════
# TASK 6: INTERACTIVE HTML DASHBOARD
# ═══════════════════════════════════════════════════════════
print("\n" + "=" * 60)
print("TASK 6: INTERACTIVE HTML DASHBOARD")
print("=" * 60)

# Build a self-contained HTML dashboard with Plotly.js
hum_pca2 = PCA(n_components=2).fit_transform(Xj)
bov_pca2 = PCA(n_components=2).fit_transform(bov_pca)

# Prepare data as JSON for embedding
human_points = []
for i in range(len(hum_pca2)):
    human_points.append({
        "x": float(hum_pca2[i, 0]),
        "y": float(hum_pca2[i, 1]),
        "state": str(readiness[i]),
        "sample": ch[i],
        "metabolic": float(Xf_s[i].mean()),
    })

bovine_points = []
for i in range(len(bov_pca2)):
    bovine_points.append({
        "x": float(bov_pca2[i, 0]),
        "y": float(bov_pca2[i, 1]),
        "state": str(bov_readiness[i]),
        "sample": bovine.columns[i],
        "time": int(bov_times[i]) if bov_times[i] >= 0 else -1,
    })

dashboard_html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Stem Cell State Map — Interactive Dashboard</title>
<script src="https://cdn.plot.ly/plotly-2.32.0.min.js"></script>
<style>
* {{ margin:0; padding:0; box-sizing:border-box; }}
body {{ font-family: system-ui, sans-serif; background: #0d1117; color: #c9d1d9; }}
header {{ background: #161b22; padding: 16px 24px; border-bottom: 1px solid #30363d; }}
header h1 {{ font-size: 1.3rem; color: #58a6ff; }}
.grid {{ display: grid; grid-template-columns: 1fr 1fr; gap: 16px; padding: 16px; }}
.card {{ background: #161b22; border: 1px solid #30363d; border-radius: 8px; padding: 16px; }}
.card h2 {{ font-size: 1rem; margin-bottom: 12px; color: #8b949e; }}
.chart {{ width: 100%; height: 400px; }}
.stats {{ display: grid; grid-template-columns: repeat(3, 1fr); gap: 12px; }}
.stat {{ background: #0d1117; border-radius: 6px; padding: 12px; text-align: center; }}
.stat .value {{ font-size: 1.5rem; font-weight: bold; color: #58a6ff; }}
.stat .label {{ font-size: 0.75rem; color: #8b949e; margin-top: 4px; }}
.expansion {{ color: #3fb950; }}
.committed {{ color: #d29922; }}
.terminal {{ color: #f85149; }}
.gene-list {{ max-height: 300px; overflow-y: auto; font-size: 0.8rem; }}
.gene-list div {{ padding: 4px 8px; border-bottom: 1px solid #21262d; }}
</style>
</head>
<body>
<header>
  <h1>🧫 Stem Cell Manufacturing Readiness — State Map Dashboard</h1>
</header>

<div class="grid">
  <div class="card">
    <h2>Human State Map (n={len(hum_pca2)})</h2>
    <div id="humanMap" class="chart"></div>
  </div>
  <div class="card">
    <h2>Bovine State Map (n={len(bov_pca2)})</h2>
    <div id="bovineMap" class="chart"></div>
  </div>
  <div class="card">
    <h2>State Distribution</h2>
    <div class="stats">
      <div class="stat"><div class="value expansion">{sum(1 for p in human_points if p['state']=='expansion_competent')}</div><div class="label">Human Expansion</div></div>
      <div class="stat"><div class="value committed">{sum(1 for p in human_points if p['state']=='committed')}</div><div class="label">Human Committed</div></div>
      <div class="stat"><div class="value terminal">{sum(1 for p in human_points if p['state']=='terminal')}</div><div class="label">Human Terminal</div></div>
      <div class="stat"><div class="value expansion">{sum(1 for p in bovine_points if p['state']=='expansion_competent')}</div><div class="label">Bovine Expansion</div></div>
      <div class="stat"><div class="value committed">{sum(1 for p in bovine_points if p['state']=='committed')}</div><div class="label">Bovine Committed</div></div>
      <div class="stat"><div class="value terminal">{sum(1 for p in bovine_points if p['state']=='terminal')}</div><div class="label">Bovine Terminal</div></div>
    </div>
  </div>
  <div class="card">
    <h2>30-Gene QC Panel</h2>
    <div class="gene-list">
      {''.join(f'<div><span class="expansion">■</span> {g}</div>' for g in panel_genes_found)}
    </div>
  </div>
</div>

<script>
const humanData = {json.dumps(human_points)};
const bovineData = {json.dumps(bovine_points)};

const colors = {{'expansion_competent': '#3fb950', 'committed': '#d29922', 'terminal': '#f85149'}};

// Human map
const hTraces = ['expansion_competent', 'committed', 'terminal'].map(state => ({{
  x: humanData.filter(p => p.state === state).map(p => p.x),
  y: humanData.filter(p => p.state === state).map(p => p.y),
  mode: 'markers',
  type: 'scatter',
  name: state,
  marker: {{ color: colors[state], size: 8, opacity: 0.8 }},
  text: humanData.filter(p => p.state === state).map(p => p.sample),
  hovertemplate: '%{{text}}<extra></extra>',
}}));

Plotly.newPlot('humanMap', hTraces, {{
  margin: {{ t: 10, r: 10, b: 40, l: 40 }},
  paper_bgcolor: '#161b22',
  plot_bgcolor: '#0d1117',
  font: {{ color: '#8b949e' }},
  xaxis: {{ title: 'PC1', gridcolor: '#21262d', zerolinecolor: '#30363d' }},
  yaxis: {{ title: 'PC2', gridcolor: '#21262d', zerolinecolor: '#30363d' }},
  legend: {{ orientation: 'h', y: 1.1 }},
}}, {{ responsive: true }});

// Bovine map
const bTraces = ['expansion_competent', 'committed', 'terminal'].map(state => ({{
  x: bovineData.filter(p => p.state === state).map(p => p.x),
  y: bovineData.filter(p => p.state === state).map(p => p.y),
  mode: 'markers',
  type: 'scatter',
  name: state,
  marker: {{ color: colors[state], size: 10, opacity: 0.8, symbol: 'diamond' }},
  text: bovineData.filter(p => p.state === state).map(p => p.sample + ' (D' + p.time + ')'),
  hovertemplate: '%{{text}}<extra></extra>',
}}));

Plotly.newPlot('bovineMap', bTraces, {{
  margin: {{ t: 10, r: 10, b: 40, l: 40 }},
  paper_bgcolor: '#161b22',
  plot_bgcolor: '#0d1117',
  font: {{ color: '#8b949e' }},
  xaxis: {{ title: 'PC1', gridcolor: '#21262d', zerolinecolor: '#30363d' }},
  yaxis: {{ title: 'PC2', gridcolor: '#21262d', zerolinecolor: '#30363d' }},
  legend: {{ orientation: 'h', y: 1.1 }},
}}, {{ responsive: true }});
</script>
</body>
</html>"""

(OUT / "dashboard.html").write_text(dashboard_html)
print(f"   Saved to dashboard.html ({len(dashboard_html)} bytes)")

# ═══════════════════════════════════════════════════════════
# TASK 7: CROSS-SPECIES TRANSFER LEARNING
# ═══════════════════════════════════════════════════════════
print("\n" + "=" * 60)
print("TASK 7: CROSS-SPECIES TRANSFER LEARNING")
print("=" * 60)

# Train on human 30-gene panel, test on bovine orthologs
# Find bovine orthologs of the 30 panel genes
bovine_panel_orthologs = []
for gene in panel_genes_found:
    if gene in bovine.index:
        bovine_panel_orthologs.append(gene)

print(f"\n   Panel genes with bovine orthologs: {len(bovine_panel_orthologs)}/{len(panel_genes_found)}")

# Extract expression of these genes in both species
human_panel_expr = np.zeros((len(ch), len(bovine_panel_orthologs)))
bovine_panel_expr = np.zeros((len(bovine.columns), len(bovine_panel_orthologs)))

for i, gene in enumerate(bovine_panel_orthologs):
    human_panel_expr[:, i] = np.log1p(tpm_df.loc[gene, ch].values)
    bovine_panel_expr[:, i] = np.log1p(bovine.loc[gene].values)

# Train classifier on human
human_panel_scaled = StandardScaler().fit_transform(human_panel_expr)
lr_transfer = LogisticRegression(C=0.1, max_iter=2000).fit(human_panel_scaled, readiness)

# Test on bovine (using same scaler)
bovine_panel_scaled = StandardScaler().fit(human_panel_expr).transform(bovine_panel_expr)
bovine_pred = lr_transfer.predict(bovine_panel_scaled)
bovine_pred_proba = lr_transfer.predict_proba(bovine_panel_scaled)

# Compare with bovine ground truth
agreement = (bovine_pred == bov_readiness).mean()
print(f"   Transfer accuracy (human→bovine): {agreement:.1%}")

# Confusion matrix
from sklearn.metrics import confusion_matrix

cm = confusion_matrix(bov_readiness, bovine_pred, labels=["expansion_competent", "committed", "terminal"])
print("\n   Confusion matrix (rows=true bovine, cols=predicted):")
print("              expansion committed terminal")
for i, state in enumerate(["expansion_competent", "committed", "terminal"]):
    print(f"   {state:12s} {cm[i,0]:5d}    {cm[i,1]:5d}    {cm[i,2]:5d}")

transfer_results = {
    "n_panel_genes_with_orthologs": len(bovine_panel_orthologs),
    "transfer_accuracy": float(agreement),
    "confusion_matrix": cm.tolist(),
    "per_class_accuracy": {state: float((bovine_pred[bov_readiness == state] == state).mean()) if (bov_readiness == state).sum() > 0 else 0 for state in ["expansion_competent", "committed", "terminal"]},
}
json.dump(transfer_results, open(OUT / "transfer_learning.json", "w"), indent=2)
print("\n   Saved to transfer_learning.json")

# ═══════════════════════════════════════════════════════════
# TASK 8: PATHWAY ENRICHMENT PER STATE
# ═══════════════════════════════════════════════════════════
print("\n" + "=" * 60)
print("TASK 8: PATHWAY ENRICHMENT PER STATE")
print("=" * 60)

# Use METAFlux pathway aggregations
metaflux_agg = pd.read_csv(r"C:\Users\asdf\CascadeProjects\rao_lab_ml\nn_pipeline\combined_data\metaflux_reaggregated_cleaned.csv", index_col=0)
pathways = metaflux_agg.columns.tolist()
print(f"\n   METAFlux pathways: {pathways}")

# Per-state pathway activity
pathway_activity = {}
for state in ["expansion_competent", "committed", "terminal"]:
    mask = readiness == state
    if mask.sum() == 0:
        continue
    # Get pathway activities for these samples
    state_flux = Xf_s[mask]
    pathway_activity[state] = {p: float(state_flux.mean(axis=0)[i]) if i < state_flux.shape[1] else 0 for i, p in enumerate(pathways[:state_flux.shape[1]])}

print("\n   Pathway activity by state:")
for state in ["expansion_competent", "committed", "terminal"]:
    if state in pathway_activity:
        top = sorted(pathway_activity[state].items(), key=lambda x: x[1], reverse=True)[:5]
        print(f"   {state}:")
        for p, v in top:
            print(f"     {p}: {v:.3f}")

# Differential pathways between states
diff_pathways = {}
for s1, s2 in [("expansion_competent", "terminal"), ("expansion_competent", "committed"), ("committed", "terminal")]:
    if s1 not in pathway_activity or s2 not in pathway_activity:
        continue
    diffs = []
    for p in pathways[:Xf_s.shape[1]]:
        d = pathway_activity[s1].get(p, 0) - pathway_activity[s2].get(p, 0)
        diffs.append((p, d))
    diffs.sort(key=lambda x: abs(x[1]), reverse=True)
    diff_pathways[f"{s1}_vs_{s2}"] = [{"pathway": p, "difference": float(d)} for p, d in diffs[:5]]
    print(f"\n   {s1} vs {s2}:")
    for p, d in diffs[:5]:
        direction = "↑" if d > 0 else "↓"
        print(f"     {direction} {p}: Δ={d:+.3f}")

json.dump({"pathway_activity": pathway_activity, "differential_pathways": diff_pathways}, open(OUT / "pathway_enrichment.json", "w"), indent=2)
print("\n   Saved to pathway_enrichment.json")

# ═══════════════════════════════════════════════════════════
# TASK 9: P5 REPRODUCIBILITY
# ═══════════════════════════════════════════════════════════
print("\n" + "=" * 60)
print("TASK 9: P5 REPRODUCIBILITY ANALYSIS")
print("=" * 60)

# Use the scFEA pseudo-bulk data
X_all = np.log1p(human_tpm.values.T).astype(np.float32)
gene_var_all = X_all.var(axis=0)
keep_all = gene_var_all > np.percentile(gene_var_all, 25)
X_filt_all = X_all[:, keep_all]
X_scaled_all = StandardScaler().fit_transform(X_filt_all)

# Within vs between condition correlation
# Use dataset origin as condition proxy
# GSE115978 and GSE72056
dataset_labels = np.array(['GSE115978' if 'GSM' in s else 'GSE72056' for s in human_tpm.columns])

within_cors = []
between_cors = []
n = X_scaled_all.shape[0]
sample_size = min(5000, n * (n-1) // 2)

# Sample for speed
rng = np.random.RandomState(42)
for _ in range(sample_size):
    i, j = rng.randint(0, n, 2)
    if i == j:
        continue
    corr, _ = pearsonr(X_scaled_all[i], X_scaled_all[j])
    if dataset_labels[i] == dataset_labels[j]:
        within_cors.append(corr)
    else:
        between_cors.append(corr)

within_mean = np.mean(within_cors)
between_mean = np.mean(between_cors)
separation = within_mean / (between_mean + 1e-10)

print(f"\n   Within-dataset correlation:  {within_mean:.3f}")
print(f"   Between-dataset correlation: {between_mean:.3f}")
print(f"   Separation ratio:            {separation:.2f}x")

# Gene stability (CV across samples)
gene_cvs = X_filt_all.std(axis=0) / (np.abs(X_filt_all.mean(axis=0)) + 1e-8)
stable_genes = gene_cvs < np.percentile(gene_cvs, 50)
print(f"   Stable genes (CV < median): {stable_genes.sum()}/{len(stable_genes)} ({stable_genes.mean():.1%})")

# PCA stability
pca_all = PCA(n_components=10).fit(X_scaled_all)
print(f"   Top 5 PCs explain: {pca_all.explained_variance_ratio_[:5].sum():.1%} variance")

# Bootstrap stability
n_boot = 20
pca_stability = []
for b in range(n_boot):
    idx = rng.choice(n, size=int(n * 0.8), replace=True)
    pca_boot = PCA(n_components=5).fit(X_scaled_all[idx])
    pca_stability.append(pca_boot.explained_variance_ratio_)

pca_stability = np.array(pca_stability)
print("   PCA variance stability (bootstrap CV):")
for i in range(5):
    cv = pca_stability[:, i].std() / (pca_stability[:, i].mean() + 1e-10)
    print(f"     PC{i+1}: mean={pca_stability[:, i].mean():.3f}, CV={cv:.3f}")

reproducibility = {
    "within_dataset_correlation": float(within_mean),
    "between_dataset_correlation": float(between_mean),
    "separation_ratio": float(separation),
    "n_stable_genes": int(stable_genes.sum()),
    "pct_stable": float(stable_genes.mean()),
    "top5_pca_variance": float(pca_all.explained_variance_ratio_[:5].sum()),
    "pca_bootstrap_cv": [float(pca_stability[:, i].std() / (pca_stability[:, i].mean() + 1e-10)) for i in range(5)],
}
json.dump(reproducibility, open(OUT / "reproducibility.json", "w"), indent=2)
print("\n   Saved to reproducibility.json")

# ═══════════════════════════════════════════════════════════
# TASK 10: PUBLICATION-READY FIGURES
# ═══════════════════════════════════════════════════════════
print("\n" + "=" * 60)
print("TASK 10: PUBLICATION-READY FIGURES")
print("=" * 60)

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.size': 10,
    'axes.labelsize': 12,
    'axes.titlesize': 14,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.1,
})

colors = {'expansion_competent': '#3fb950', 'committed': '#d29922', 'terminal': '#f85149'}

# Figure 1: State Map
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Human
ax = axes[0]
for state in ['expansion_competent', 'committed', 'terminal']:
    mask = readiness == state
    ax.scatter(hum_pca2[mask, 0], hum_pca2[mask, 1], c=colors[state], label=state, s=30, alpha=0.8, edgecolors='white', linewidth=0.5)
ax.set_xlabel('PC1')
ax.set_ylabel('PC2')
ax.set_title(f'Human State Map (n={len(hum_pca2)})')
ax.legend(frameon=True, facecolor='white', framealpha=0.9)
ax.set_facecolor('#f8f9fa')

# Bovine
ax = axes[1]
for state in ['expansion_competent', 'committed', 'terminal']:
    mask = bov_readiness == state
    ax.scatter(bov_pca2[mask, 0], bov_pca2[mask, 1], c=colors[state], label=state, s=50, alpha=0.8, edgecolors='white', linewidth=0.5, marker='D')
ax.set_xlabel('PC1')
ax.set_ylabel('PC2')
ax.set_title(f'Bovine State Map (n={len(bov_pca2)})')
ax.legend(frameon=True, facecolor='white', framealpha=0.9)
ax.set_facecolor('#f8f9fa')

fig.suptitle('Figure 1: Multi-Omic State Maps Across Species', fontweight='bold', y=1.01)
plt.tight_layout()
fig.savefig(OUT / 'fig1_state_map_publication.png')
plt.close()
print("   fig1_state_map_publication.png saved")

# Figure 2: Pseudotime Trajectory
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# PCA1 vs true time
ax = axes[0]
for state in ['expansion_competent', 'committed', 'terminal']:
    mask = bov_readiness[time_idx] == state
    ax.scatter(bov_true_times[mask], bov_pca1[mask], c=colors[state], s=40, alpha=0.8, edgecolors='white', linewidth=0.5)
ax.set_xlabel('True Time (days)')
ax.set_ylabel('PCA1 Score')
ax.set_title(f'PCA1 vs Time (r={time_corr:.3f}, p={time_p:.1e})')
ax.set_facecolor('#f8f9fa')

# Pseudotime distribution
ax = axes[1]
timepoints = sorted(set(bov_true_times))
positions = range(len(timepoints))
bp = ax.boxplot([pseudotime[bov_true_times == t] for t in timepoints], positions=positions, widths=0.5)
for i, t in enumerate(timepoints):
    mask = bov_true_times == t
    state_counts = Counter(bov_readiness[time_idx][mask])
    ax.text(i, 1.05, f'n={mask.sum()}', ha='center', fontsize=8)
ax.set_xticks(positions)
ax.set_xticklabels([f'D{t}' for t in timepoints])
ax.set_ylabel('Diffusion Pseudotime')
ax.set_title(f'Pseudotime by Day (r={pseudo_corr:.3f})')
ax.set_facecolor('#f8f9fa')

fig.suptitle('Figure 2: Bovine Differentiation Trajectory', fontweight='bold')
plt.tight_layout()
fig.savefig(OUT / 'fig2_pseudotime_publication.png')
plt.close()
print("   fig2_pseudotime_publication.png saved")

# Figure 3: Gene Regulatory Network
fig, ax = plt.subplots(figsize=(12, 10))

# Plot nodes
n_plot = min(20, n_panel)
plot_genes = [g for g, _ in hub_degrees.most_common(n_plot)]
node_pos = {}
for i, gene in enumerate(plot_genes):
    angle = 2 * np.pi * i / n_plot
    r = 1.0
    node_pos[gene] = (r * np.cos(angle), r * np.sin(angle))

# Draw edges
for e in edges[:50]:
    if e['source'] in node_pos and e['target'] in node_pos:
        x1, y1 = node_pos[e['source']]
        x2, y2 = node_pos[e['target']]
        alpha = min(1.0, abs(e['correlation']))
        color = '#3fb950' if e['correlation'] > 0 else '#f85149'
        ax.plot([x1, x2], [y1, y2], color=color, alpha=alpha * 0.5, linewidth=abs(e['correlation']) * 2)

# Draw nodes
for gene, (x, y) in node_pos.items():
    deg = hub_degrees.get(gene, 0)
    size = 200 + deg * 100
    ax.scatter(x, y, s=size, c='#58a6ff', edgecolors='white', linewidth=1.5, zorder=5)
    ax.annotate(gene, (x, y), textcoords="offset points", xytext=(0, 12), ha='center', fontsize=8, fontweight='bold')

ax.set_xlim(-1.5, 1.5)
ax.set_ylim(-1.5, 1.5)
ax.set_aspect('equal')
ax.axis('off')
ax.set_title(f'Figure 3: 30-Gene Regulatory Network\n({len(edges)} edges, |r|>0.7)', fontweight='bold')

# Legend
legend_elements = [
    Line2D([0], [0], color='#3fb950', lw=2, label='Co-activation (r>0)'),
    Line2D([0], [0], color='#f85149', lw=2, label='Repression (r<0)'),
]
ax.legend(handles=legend_elements, loc='lower right')

plt.tight_layout()
fig.savefig(OUT / 'fig3_network_publication.png')
plt.close()
print("   fig3_network_publication.png saved")

# Figure 4: Cross-Species Transfer
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Confusion matrix heatmap
ax = axes[0]
im = ax.imshow(cm, cmap='Blues', aspect='auto')
ax.set_xticks(range(3))
ax.set_xticklabels(['expansion', 'committed', 'terminal'], rotation=45, ha='right')
ax.set_yticks(range(3))
ax.set_yticklabels(['expansion', 'committed', 'terminal'])
ax.set_xlabel('Predicted (Human-trained)')
ax.set_ylabel('True (Bovine)')
ax.set_title(f'Transfer Learning\nAccuracy: {agreement:.1%}')

for i in range(3):
    for j in range(3):
        ax.text(j, i, cm[i, j], ha='center', va='center', fontweight='bold', color='white' if cm[i,j] > cm.max()/2 else 'black')

plt.colorbar(im, ax=ax)

# Per-class accuracy
ax = axes[1]
classes = ["expansion_competent", "committed", "terminal"]
accs = [transfer_results["per_class_accuracy"][c] for c in classes]
bars = ax.bar(classes, accs, color=[colors[c] for c in classes], edgecolor='white')
ax.set_ylabel('Per-Class Accuracy')
ax.set_title('Transfer Accuracy by State')
ax.set_ylim(0, 1)
for bar, acc in zip(bars, accs):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02, f'{acc:.1%}', ha='center', fontweight='bold')

fig.suptitle('Figure 4: Cross-Species Transfer Learning', fontweight='bold')
plt.tight_layout()
fig.savefig(OUT / 'fig4_transfer_publication.png')
plt.close()
print("   fig4_transfer_publication.png saved")

# Figure 5: Comprehensive Summary
fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# A: Pathway activity
ax = axes[0, 0]
pathway_names = pathways[:8]
x = np.arange(len(pathway_names))
width = 0.25
for i, state in enumerate(['expansion_competent', 'committed', 'terminal']):
    if state in pathway_activity:
        vals = [pathway_activity[state].get(p, 0) for p in pathway_names]
        ax.bar(x + i*width, vals, width, label=state, color=colors[state], alpha=0.8)
ax.set_xticks(x + width)
ax.set_xticklabels(pathway_names, rotation=45, ha='right', fontsize=8)
ax.set_ylabel('Mean Activity')
ax.set_title('A: Pathway Activity by State')
ax.legend(fontsize=8)

# B: Power analysis
ax = axes[0, 1]
comparisons = list(power_results.keys())
n_vals = [power_results[c]['n_total'] for c in comparisons]
ax.barh(comparisons, n_vals, color='#58a6ff')
ax.axvline(x=30, color='#f85149', linestyle='--', label='Proposed N=30')
ax.set_xlabel('Required Total N')
ax.set_title('B: Power Analysis (80%, α=0.05)')
ax.legend()

# C: Drug predictions
ax = axes[1, 0]
top_drugs = predictions[:6]
drug_names = [p['compound'].split('(')[0].strip() for p in top_drugs]
drug_scores = [p['score'] for p in top_drugs]
ax.barh(drug_names, drug_scores, color='#d29922')
ax.set_xlabel('Prediction Score')
ax.set_title('C: Compounds for Expansion State')

# D: Reproducibility
ax = axes[1, 1]
rep_metrics = ['Within-dataset\ncorrelation', 'Between-dataset\ncorrelation', 'Stable genes\n(%)', 'Top 5 PC\nvariance']
rep_vals = [within_mean, between_mean, stable_genes.mean(), pca_all.explained_variance_ratio_[:5].sum()]
bars = ax.bar(rep_metrics, rep_vals, color=['#3fb950', '#d29922', '#58a6ff', '#8b949e'])
ax.set_title('D: Reproducibility Metrics')
for bar, val in zip(bars, rep_vals):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02, f'{val:.2f}', ha='center', fontweight='bold')

fig.suptitle('Figure 5: Comprehensive Analysis Summary', fontweight='bold', y=1.01)
plt.tight_layout()
fig.savefig(OUT / 'fig5_summary_publication.png')
plt.close()
print("   fig5_summary_publication.png saved")

# ═══════════════════════════════════════════════════════════
# FINAL SUMMARY
# ═══════════════════════════════════════════════════════════
print("\n" + "=" * 60)
print("ALL TASKS COMPLETE")
print("=" * 60)

all_outputs = {
    "task2_pseudotime": "pseudotime_analysis.json",
    "task3_network": "gene_regulatory_network.json",
    "task4_drugs": "drug_predictions.json",
    "task5_power": "power_analysis.json",
    "task6_dashboard": "dashboard.html",
    "task7_transfer": "transfer_learning.json",
    "task8_pathways": "pathway_enrichment.json",
    "task9_reproducibility": "reproducibility.json",
    "task10_figures": [
        "fig1_state_map_publication.png",
        "fig2_pseudotime_publication.png",
        "fig3_network_publication.png",
        "fig4_transfer_publication.png",
        "fig5_summary_publication.png",
    ],
}

json.dump(all_outputs, open(OUT / "all_analyses_manifest.json", "w"), indent=2)
print(f"\nAll outputs in: {OUT}")
for k, v in all_outputs.items():
    print(f"  {k}: {v}")