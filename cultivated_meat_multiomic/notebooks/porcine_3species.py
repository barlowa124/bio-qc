"""Porcine GSE206914 processing + 3-species comparison."""
import json
import re
import time
import urllib.request
import warnings
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")
PROJ = Path(__file__).resolve().parents[1]
CROSS = PROJ / "cross_species_validation"
OUT = PROJ / "p2_state_map/output"
OUT.mkdir(exist_ok=True)

print("=" * 60)
print("PORCINE GSE206914 + 3-SPECIES COMPARISON")
print("=" * 60)

# Load human data
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
readiness = np.array([rmap[c] for c in cl])

print(f"Human: {len(ch)} samples, states={dict(Counter(readiness))}")

# ── Check porcine datasets ──
porcine_info = []
for gse_id in ['GSE206912', 'GSE206913', 'GSE206914']:
    print(f"\nChecking {gse_id}...")
    try:
        req = urllib.request.Request(
            f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse_id}",
            headers={"User-Agent": "RaoLab/1.0"}
        )
        with urllib.request.urlopen(req, timeout=30) as r:
            text = r.read().decode("utf-8", errors="ignore")
        title_m = re.search(r'<td[^>]*>Title</td>.*?<td[^>]*>(.*?)</td>', text, re.DOTALL)
        title = re.sub(r'<[^>]+>', '', title_m.group(1)).strip() if title_m else "N/A"
        org_m = re.search(r'<td[^>]*>Organism</td>.*?<td[^>]*>(.*?)</td>', text, re.DOTALL)
        org = re.sub(r'<[^>]+>', '', org_m.group(1)).strip() if org_m else "N/A"
        gsms = list(set(re.findall(r'GSM\d+', text)))
        suppl = re.findall(r'https?://[^\s\"\'<>]+\.(?:txt|gz|csv|tsv|zip|mtx|h5|h5ad)[^\s\"\'<>]*', text)
        info = {"gse": gse_id, "title": title[:200], "organism": org, "n_samples": len(gsms), "n_suppl": len(suppl)}
        porcine_info.append(info)
        print(f"  {title[:120]}")
        print(f"  Organism={org}, Samples={len(gsms)}, Suppl={len(suppl)}")
        for u in suppl[:3]:
            print(f"  Suppl: {u[-100:]}")
    except Exception as e:
        print(f"  Error: {e}")
        porcine_info.append({"gse": gse_id, "error": str(e)})
    time.sleep(1.0)

# ── 3-species marker comparison ──
print("\n─── 3-Species Marker Comparison ───")

conserved = ["PAX7","MYF5","MYOD1","MYOG","MYH3","DES","MYL1","NEB",
             "POU5F1","SOX2","NANOG","MYC","KLF4","CD44","CD90","GAPDH","ACTB"]

# Human markers
hum_log = np.log1p(tpm_df[ch].values.T).astype(np.float32)
human_markers = {}
for gene in conserved:
    if gene in tpm_df.index:
        gi = list(tpm_df.index).index(gene)
        expr = hum_log[:, gi]
        human_markers[gene] = {
            "mean": float(expr.mean()),
            "by_state": {s: float(expr[readiness == s].mean()) for s in ["expansion_competent","committed","terminal"]}
        }

# Bovine markers from existing results
bov_data = json.loads((OUT / "cross_species_comparison.json").read_text())
bovine_markers = bov_data["cross_species_markers"]["bovine"]

# Porcine ortholog mapping
pig_orthologs = {g: g for g in conserved}  # Most myogenic genes have identical symbols

# Build comparison table
comparison_rows = []
for gene in conserved:
    row = {"gene": gene}
    if gene in human_markers:
        row["human"] = human_markers[gene]["by_state"]
        row["human_mean"] = human_markers[gene]["mean"]
    if gene in bovine_markers:
        row["bovine"] = bovine_markers[gene]["by_state"]
        row["bovine_mean"] = bovine_markers[gene]["mean"]
    row["porcine_ortholog"] = pig_orthologs.get(gene, "?")
    comparison_rows.append(row)

# Compute cross-species correlations
human_vecs = []
bovine_vecs = []
for row in comparison_rows:
    if "human" in row and "bovine" in row:
        h = row["human"]
        b = row["bovine"]
        for s in ["expansion_competent","committed","terminal"]:
            human_vecs.append(h[s])
            bovine_vecs.append(b[s])

if len(human_vecs) > 2:
    r, p = pearsonr(human_vecs, bovine_vecs)
    print(f"\nHuman-Bovine marker correlation: r={r:.3f}, p={p:.4f}")

# Print table
print(f"\n{'Gene':12s} {'Human (E/C/T)':30s} {'Bovine (E/C/T)':30s} {'Porcine':10s}")
print("-" * 85)
for row in comparison_rows:
    h_str = f"{row['human']['expansion_competent']:.2f}/{row['human']['committed']:.2f}/{row['human']['terminal']:.2f}" if "human" in row else "N/A"
    b_str = f"{row['bovine']['expansion_competent']:.2f}/{row['bovine']['committed']:.2f}/{row['bovine']['terminal']:.2f}" if "bovine" in row else "N/A"
    p_str = row.get("porcine_ortholog", "?")
    print(f"{row['gene']:12s} {h_str:30s} {b_str:30s} {p_str:10s}")

# ── Save ──
results = {
    "porcine_datasets": porcine_info,
    "three_species_comparison": comparison_rows,
    "human_bovine_correlation": float(r) if len(human_vecs) > 2 else None,
    "n_conserved_markers": len([r for r in comparison_rows if "human" in r and "bovine" in r]),
    "porcine_status": "metadata_only" if not any("count" in str(p) for p in porcine_info) else "data_available"
}
json.dump(results, open(OUT / "three_species_comparison.json", "w"), indent=2)
print("\nSaved to three_species_comparison.json")
print("=" * 60)
print("DONE")
