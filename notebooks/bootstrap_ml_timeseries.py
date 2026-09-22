"""Bootstrap stability, ML comparison, time-series modeling."""
import json
import warnings
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC

warnings.filterwarnings("ignore")
PROJ = Path(__file__).resolve().parents[1]
CROSS = PROJ / "cross_species_validation"
OUT = PROJ / "p2_state_map/output"
OUT.mkdir(exist_ok=True)

print("=" * 60)
print("BOOTSTRAP + ML COMPARISON + TIME-SERIES")
print("=" * 60)

# Load data
tpm_df = pd.read_csv(r"C:\Users\asdf\CascadeProjects\rao_lab_ml\nn_pipeline\combined_data\metaflux_scfea_pseudobulk\scfea_pseudobulk_tpm.csv", index_col=0)
flux_df = pd.read_csv(r"C:\Users\asdf\CascadeProjects\rao_lab_ml\nn_pipeline\combined_data\metaflux_scfea_pseudobulk\metaflux_combined_flux_matrix_all.csv", index_col=0)
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
y = np.array([rmap[c] for c in cl])

# 30-gene panel
qc_panel_path = PROJ / "p3_qc_panel/output/qc_panel_239.json"
if qc_panel_path.exists():
    qc_data = json.loads(qc_panel_path.read_text())
    gp = qc_data.get("panel_genes", qc_data.get("genes", []))
    if isinstance(gp[0], dict):
        gp = [g["gene"] if isinstance(g, dict) else g for g in gp]
else:
    gp = ["UXS1","PLOD1","MALAT1","C1D","KIF1B","COL1A1","FN1","VIM","ACTA2","MYH9",
          "TPM1","TAGLN","CNN1","MYL9","ACTG2","DES","MYH11","LMOD1","PDLIM3","MYLK",
          "ITGA8","SYNPO2","MRVI1","PTGIS","GUCY1A1","NEXN","PPP1R12B","SORBS1","FLNC","SPARCL1"]

# Build gene expression matrix for panel genes
panel_genes_in_data = [g for g in gp if g in tpm_df.index]
print(f"Panel genes in data: {len(panel_genes_in_data)}/{len(gp)}")

X_panel = np.zeros((len(ch), len(panel_genes_in_data)), dtype=np.float32)
for i, g in enumerate(panel_genes_in_data):
    gi = list(tpm_df.index).index(g)
    X_panel[:, i] = Xr[:, gi]

# ── Bootstrap stability ──
print("\n─── Bootstrap Stability (1000 iterations) ───")
n_boot = 1000
boot_scores = []
boot_genes = []
rng = np.random.RandomState(42)

for b in range(n_boot):
    idx = rng.choice(len(ch), size=len(ch), replace=True)
    Xb = X_panel[idx]
    yb = y[idx]
    # Only run if all 3 classes present
    if len(set(yb)) < 3:
        continue
    lr = LogisticRegression(C=1.0, penalty='l1', solver='saga', max_iter=2000, random_state=b)
    lr.fit(Xb, yb)
    s = lr.score(X_panel, y)
    boot_scores.append(s)
    coef = np.abs(lr.coef_).sum(axis=0)
    top10 = np.argsort(coef)[-10:][::-1]
    boot_genes.append([panel_genes_in_data[i] for i in top10])

boot_scores = np.array(boot_scores)
print(f"Bootstrap accuracy: {boot_scores.mean():.3f} +/- {boot_scores.std():.3f}")
print(f"95% CI: [{np.percentile(boot_scores, 2.5):.3f}, {np.percentile(boot_scores, 97.5):.3f}]")

# Gene stability
gene_freq = Counter()
for bg in boot_genes:
    for g in bg:
        gene_freq[g] += 1

print(f"\nGene selection stability (top-10 frequency over {len(boot_genes)} bootstraps):")
for gene, freq in gene_freq.most_common(15):
    print(f"  {gene}: {freq}/{len(boot_genes)} ({freq/len(boot_genes):.1%})")

# ── ML Model Comparison ──
print("\n─── ML Model Comparison ───")
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

models = {
    "LogisticRegression": LogisticRegression(C=1.0, max_iter=2000, random_state=42),
    "RandomForest": RandomForestClassifier(n_estimators=100, random_state=42),
    "SVM_RBF": SVC(kernel='rbf', random_state=42),
    "SVM_Linear": SVC(kernel='linear', random_state=42),
}

ml_results = {}
for name, model in models.items():
    scores = cross_val_score(model, X_panel, y, cv=cv, scoring='accuracy')
    ml_results[name] = {"mean": float(scores.mean()), "std": float(scores.std()), "scores": [float(s) for s in scores]}
    print(f"  {name:20s}: {scores.mean():.3f} +/- {scores.std():.3f}")

# ── Time-series modeling of bovine D0-D7 ──
print("\n─── Bovine D0-D7 Time-Series Modeling ──")

# Load bovine timecourse
tc = pd.read_csv(CROSS / "GSE173199/timecourse_counts.csv", index_col=0)
symbol_map = json.loads((CROSS / "bovine_ensembl_to_symbol.json").read_text())
new_index = [symbol_map.get(g.split(".")[0], g) for g in tc.index.tolist()]
tc.index = new_index
tc = tc[~tc.index.duplicated(keep='first')]

# Parse timepoints from column names
timepoints = []
for col in tc.columns:
    # Extract day number
    import re
    m = re.search(r'[Dd](\d+)', col)
    if m:
        timepoints.append(int(m.group(1)))
    else:
        timepoints.append(-1)

tc_times = np.array(timepoints)
valid = tc_times >= 0
tc_valid = tc.loc[:, [tc.columns[i] for i in range(len(tc.columns)) if valid[i]]]
tc_times = tc_times[valid]

print(f"Timepoints: {sorted(set(tc_times))}")
print(f"Samples: {len(tc_times)}")

# Normalize
tc_cpm = tc_valid.values / tc_valid.values.sum(axis=0) * 1e6
tc_log = np.log1p(tc_cpm.T).astype(np.float32)

# Temporal clustering: group samples by timepoint
time_groups = {}
for i, t in enumerate(tc_times):
    if t not in time_groups:
        time_groups[t] = []
    time_groups[t].append(i)

# Compute mean expression per timepoint
time_means = {}
for t, idxs in sorted(time_groups.items()):
    time_means[t] = tc_log[idxs].mean(axis=0)

# Find genes with monotonic trends
trend_scores = {}
for gi, gene in enumerate(tc_valid.index):
    expr = [time_means[t][gi] for t in sorted(time_means.keys())]
    times = sorted(time_means.keys())
    if len(expr) > 2:
        r, p = pearsonr(times, expr)
        trend_scores[gene] = {"r": float(r), "p": float(p), "expression": [float(e) for e in expr]}

# Top trending genes
trending = sorted(trend_scores.items(), key=lambda x: abs(x[1]["r"]), reverse=True)
print("\nTop temporally-correlated genes:")
for gene, info in trending[:10]:
    print(f"  {gene}: r={info['r']:.3f}, p={info['p']:.4f}")

# Check panel gene trends in bovine
print("\nPanel gene temporal trends in bovine:")
panel_trends = {}
for gene in gp:
    if gene in trend_scores:
        panel_trends[gene] = trend_scores[gene]
        print(f"  {gene}: r={trend_scores[gene]['r']:.3f}, p={trend_scores[gene]['p']:.4f}")

# ── Save ──
def convert(obj):
    if isinstance(obj, (np.integer,)): return int(obj)
    if isinstance(obj, (np.floating,)): return float(obj) if not np.isnan(obj) else None
    if isinstance(obj, np.ndarray): return obj.tolist()
    return obj

results = {
    "bootstrap": {
        "n_iterations": int(len(boot_scores)),
        "mean_accuracy": float(boot_scores.mean()),
        "std_accuracy": float(boot_scores.std()),
        "ci_95": [float(np.percentile(boot_scores, 2.5)), float(np.percentile(boot_scores, 97.5))],
        "gene_stability": {g: float(f)/len(boot_genes) for g, f in gene_freq.most_common(30)}
    },
    "ml_comparison": ml_results,
    "time_series": {
        "n_timepoints": int(len(time_means)),
        "timepoints": [int(t) for t in sorted(time_means.keys())],
        "n_genes_tested": int(len(trend_scores)),
        "top_trending": [{"gene": str(g), "r": float(i["r"]) if not np.isnan(i["r"]) else None, "p": float(i["p"]) if not np.isnan(i["p"]) else None} for g, i in trending[:20]],
        "panel_gene_trends": {str(g): {"r": float(i["r"]) if not np.isnan(i["r"]) else None, "p": float(i["p"]) if not np.isnan(i["p"]) else None} for g, i in panel_trends.items()}
    }
}
json.dump(results, open(OUT / "bootstrap_ml_timeseries.json", "w"), indent=2, default=convert)
print("\nSaved to bootstrap_ml_timeseries.json")
print("=" * 60)
print("DONE")
