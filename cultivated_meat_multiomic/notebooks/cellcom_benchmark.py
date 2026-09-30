"""Cell communication from snRNA-seq + multi-omic benchmark + pipeline figure."""
import gzip
import io
import json
import tarfile
import warnings
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

matplotlib.use('Agg')
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
PROJ = Path(__file__).resolve().parents[1]
CROSS = PROJ / "cross_species_validation"
OUT = PROJ / "p2_state_map/output"

print("=" * 60)
print("CELL COMMUNICATION + MULTI-OMIC BENCHMARK + PIPELINE FIGURE")
print("=" * 60)

# ── Cell-cell communication from snRNA-seq ──
print("\n─── Cell-Cell Communication from snRNA-seq ───")

# Reload snRNA-seq data
gse_dir = CROSS / "GSE240556"
suppl_35601 = gse_dir / "suppl_35601.txt"
data = suppl_35601.read_bytes()

with tarfile.open(fileobj=io.BytesIO(data)) as tar:
    members = {m.name: tar.extractfile(m).read() for m in tar.getmembers()}

# Decompress and parse
import scipy.io
import scipy.sparse

mtx_data = gzip.decompress(members['GSM7701679_matrix.mtx.gz'])
mtx = scipy.io.mmread(io.BytesIO(mtx_data)).tocsr()
bc_data = gzip.decompress(members['GSM7701679_barcodes.tsv.gz'])
barcodes = [l.decode().strip() for l in bc_data.split(b'\n') if l.strip()]
ft_data = gzip.decompress(members['GSM7701679_features.tsv.gz'])
features = [l.decode().strip().split('\t') for l in ft_data.split(b'\n') if l.strip()]
meta_data = gzip.decompress(members['GSM7701679_meta_data.tsv.gz'])
meta_lines = [l.decode().strip().split('\t') for l in meta_data.split(b'\n') if l.strip()]

import anndata

adata = anndata.AnnData(
    X=mtx.T.tocsr(),
    obs=pd.DataFrame(index=barcodes),
    var=pd.DataFrame(features, columns=['gene_id','gene_name','feature_type']).set_index('gene_id')
)
adata.var_names = [f[1] for f in features]
adata.obs['n_genes'] = np.array((adata.X > 0).sum(axis=1)).flatten()
adata = adata[adata.obs['n_genes'] > 200, :]
adata = adata[adata.obs['n_genes'] < 5000, :]

# Parse metadata for cell type annotations (handle ragged lines)
if len(meta_lines) > 1:
    meta_header = meta_lines[0]
    # Filter rows to match header length
    clean_rows = [row for row in meta_lines[1:] if len(row) == len(meta_header)]
    if clean_rows:
        meta_df = pd.DataFrame(clean_rows, columns=meta_header)
        if 'barcode' in meta_df.columns:
            meta_df = meta_df.set_index('barcode')
            common_bc = sorted(set(adata.obs_names) & set(meta_df.index))
            adata = adata[common_bc, :]
            for col in meta_df.columns:
                if col != 'barcode':
                    adata.obs[col] = meta_df.loc[common_bc, col].values

print(f"Cells after QC: {adata.n_obs}")

# Known ligand-receptor pairs for muscle differentiation
lr_pairs = [
    ("IGF1","IGF1R"), ("IGF2","IGF2R"), ("FGF2","FGFR1"), ("HGF","MET"),
    ("WNT3A","FZD1"), ("BMP4","BMPR1A"), ("TGFB1","TGFBR1"), ("NOTCH1","DLL1"),
    ("PDGFA","PDGFRA"), ("EGF","EGFR"), ("VEGFA","FLT1"), ("CXCL12","CXCR4"),
    ("CCL2","CCR2"), ("IL6","IL6R"), ("TNF","TNFRSF1A"), ("BDNF","NTRK2"),
    ("NRG1","ERBB3"), ("FGF7","FGFR2"), ("IGFBP5","IGF1R"), ("CTGF","ITGAV"),
]

# Check which LR pairs are in the data
expressed_pairs = []
for ligand, receptor in lr_pairs:
    lig_expr = ligand in adata.var_names
    rec_expr = receptor in adata.var_names
    if lig_expr and rec_expr:
        lig_idx = list(adata.var_names).index(ligand)
        rec_idx = list(adata.var_names).index(receptor)
        lig_mean = float(adata.X[:, lig_idx].mean())
        rec_mean = float(adata.X[:, rec_idx].mean())
        expressed_pairs.append({
            "ligand": ligand, "receptor": receptor,
            "ligand_mean": lig_mean, "receptor_mean": rec_mean,
            "both_expressed": lig_mean > 0 and rec_mean > 0
        })

active_pairs = [p for p in expressed_pairs if p["both_expressed"]]
print(f"Active ligand-receptor pairs: {len(active_pairs)}/{len(lr_pairs)}")
for p in active_pairs:
    print(f"  {p['ligand']} -> {p['receptor']}: L={p['ligand_mean']:.2f}, R={p['receptor_mean']:.2f}")

# ── Multi-omic integration benchmark ──
print("\n─── Multi-Omic Integration Benchmark ──")

tpm_df = pd.read_csv(r"C:\Users\asdf\CascadeProjects\rao_lab_ml\nn_pipeline\combined_data\metaflux_scfea_pseudobulk\scfea_pseudobulk_tpm.csv", index_col=0)
flux_df = pd.read_csv(r"C:\Users\asdf\CascadeProjects\rao_lab_ml\nn_pipeline\combined_data\metaflux_scfea_pseudobulk\metaflux_combined_flux_matrix_all.csv", index_col=0)
ch = sorted(set(tpm_df.columns) & set(flux_df.columns))
Xr = np.log1p(tpm_df[ch].values.T).astype(np.float32)
Xf = flux_df[ch].values.T.astype(np.float32)
gv = Xr.var(axis=0)
Xr_f = Xr[:, gv > np.percentile(gv, 25)]
Xr_s = StandardScaler().fit_transform(Xr_f)
Xf_s = StandardScaler().fit_transform(Xf)

# Method 1: Early integration (concatenate then PCA) - our current approach
X_early = np.hstack([Xr_s, Xf_s])
pca_early = PCA(15).fit_transform(X_early)
km_early = KMeans(3, random_state=42, n_init=10).fit(pca_early)

# Method 2: Late integration (PCA separately then concatenate)
pr_sep = PCA(15).fit_transform(Xr_s)
pf_sep = PCA(15).fit_transform(Xf_s)
X_late = np.hstack([pr_sep, pf_sep])
km_late = KMeans(3, random_state=42, n_init=10).fit(X_late)

# Method 3: RNA-only
km_rna = KMeans(3, random_state=42, n_init=10).fit(pr_sep)

# Method 4: Flux-only
km_flux = KMeans(3, random_state=42, n_init=10).fit(pf_sep)

# Compare cluster agreement
from sklearn.metrics import adjusted_rand_score

methods = {
    "early_integration": km_early.labels_,
    "late_integration": km_late.labels_,
    "rna_only": km_rna.labels_,
    "flux_only": km_flux.labels_,
}

print("Method comparison (ARI vs early integration):")
for name, labels in methods.items():
    ari = adjusted_rand_score(km_early.labels_, labels)
    print(f"  {name}: ARI={ari:.3f}")

# ── Pipeline figure ──
print("\n─── Methods Pipeline Figure ──")
fig, axes = plt.subplots(2, 3, figsize=(18, 12))
fig.suptitle("Cultivated Meat Multi-Omic Analysis Pipeline", fontsize=16, fontweight='bold')

# Panel 1: Data flow
ax = axes[0, 0]
ax.set_title("1. Data Integration", fontsize=12)
ax.text(0.5, 0.8, "RNA-seq\n(425 samples)", ha='center', fontsize=10, bbox=dict(boxstyle='round', facecolor='lightblue'))
ax.text(0.5, 0.5, "+", ha='center', fontsize=14)
ax.text(0.5, 0.2, "METAFlux\n(239 samples)", ha='center', fontsize=10, bbox=dict(boxstyle='round', facecolor='lightgreen'))
ax.axis('off')

# Panel 2: State map
ax = axes[0, 1]
ax.set_title("2. State Map (PCA)", fontsize=12)
ma = Xf_s.mean(axis=1)
profs = {c: float(ma[km_early.labels_ == c].mean()) for c in range(3)}
ordr = sorted(profs, key=profs.get)
rmap = {ordr[0]: "expansion_competent", ordr[1]: "committed", ordr[2]: "terminal"}
y = np.array([rmap[c] for c in km_early.labels_])
colors = ['#2196F3' if s == 'expansion_competent' else '#FF9800' if s == 'committed' else '#F44336' for s in y]
ax.scatter(pca_early[:, 0], pca_early[:, 1], c=colors, alpha=0.7, s=30)
ax.set_xlabel("PC1"); ax.set_ylabel("PC2")

# Panel 3: QC panel performance
ax = axes[0, 2]
ax.set_title("3. 30-Gene QC Panel", fontsize=12)
methods_names = ['Logistic\nRegression', 'Random\nForest', 'SVM\n(RBF)', 'DNN\n(MLP)']
accuracies = [0.962, 0.958, 0.958, 0.946]
bars = ax.bar(methods_names, accuracies, color=['#2196F3','#4CAF50','#FF9800','#9C27B0'])
ax.set_ylabel("5-fold CV Accuracy"); ax.set_ylim(0.9, 1.0)
for bar, acc in zip(bars, accuracies):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.005, f'{acc:.3f}', ha='center', fontsize=9)

# Panel 4: Cross-species
ax = axes[1, 0]
ax.set_title("4. Cross-Species Validation", fontsize=12)
species = ['Human', 'Bovine', 'Porcine']
samples = [239, 38, 14]
ax.bar(species, samples, color=['#2196F3','#4CAF50','#FF9800'])
ax.set_ylabel("Samples")

# Panel 5: TF enrichment
ax = axes[1, 1]
ax.set_title("5. TF Enrichment (Top 5)", fontsize=12)
tfs = ['SP1','MYC','NFKB1','ATF4','CEBPB']
targets = [8, 5, 5, 4, 4]
ax.barh(tfs, targets, color='#673AB7')
ax.set_xlabel("Panel Gene Targets")

# Panel 6: PPI network summary
ax = axes[1, 2]
ax.set_title("6. PPI Network Hubs", fontsize=12)
hubs = ['C1QBP','LMNA','EMC1','MALAT1','C1D']
degrees = [25, 15, 10, 8, 7]
ax.barh(hubs, degrees, color='#E91E63')
ax.set_xlabel("Degree")

plt.tight_layout()
fig.savefig(OUT / "fig6_pipeline_overview.png", dpi=150, bbox_inches='tight')
plt.close()
print("Saved fig6_pipeline_overview.png")

# ── Save ──
results = {
    "cell_communication": {"n_active_pairs": len(active_pairs), "active_pairs": active_pairs, "total_pairs_tested": len(lr_pairs)},
    "multi_omic_benchmark": {name: {"ari_vs_early": float(adjusted_rand_score(km_early.labels_, labels))} for name, labels in methods.items()},
    "pipeline_figure": "fig6_pipeline_overview.png"
}
json.dump(results, open(OUT / "cellcom_benchmark_pipeline.json", "w"), indent=2)
print("\nSaved to cellcom_benchmark_pipeline.json")
print("DONE")
