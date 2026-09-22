"""Literature benchmark + drug prediction + minimal panel + noise simulation."""
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.feature_selection import SequentialFeatureSelector
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_score
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")
PROJ = Path(__file__).resolve().parents[1]
OUT = PROJ / "p2_state_map/output"

print("=" * 70)
print("LITERATURE BENCHMARK + DRUG PREDICTION + MINIMAL PANEL + NOISE SIM")
print("=" * 70)

# ── Shared data loading ──
tpm_df = pd.read_csv(r"C:\Users\asdf\CascadeProjects\rao_lab_ml\nn_pipeline\combined_data\metaflux_scfea_pseudobulk\scfea_pseudobulk_tpm.csv", index_col=0)
flux_df = pd.read_csv(r"C:\Users\asdf\CascadeProjects\rao_lab_ml\nn_pipeline\combined_data\metaflux_scfea_pseudobulk\metaflux_combined_flux_matrix_all.csv", index_col=0)
ch = sorted(set(tpm_df.columns) & set(flux_df.columns))
Xr = np.log1p(tpm_df[ch].values.T).astype(np.float32)
Xf = flux_df[ch].values.T.astype(np.float32)
gv = Xr.var(axis=0)
Xr_f = Xr[:, gv > np.percentile(gv, 25)]
Xr_s = StandardScaler().fit_transform(Xr_f)
Xf_s = StandardScaler().fit_transform(Xf)
pr = PCA(15).fit_transform(Xr_s)
pf = PCA(15).fit_transform(Xf_s)
km = KMeans(3, random_state=42, n_init=10).fit(np.hstack([pr, pf]))
ma = Xf_s.mean(axis=1)
profs = {c: float(ma[km.labels_ == c].mean()) for c in range(3)}
ordr = sorted(profs, key=profs.get)
rmap = {ordr[0]: "expansion_competent", ordr[1]: "committed", ordr[2]: "terminal"}
y = np.array([rmap[c] for c in km.labels_])

qc_data = json.loads((PROJ / "p3_qc_panel/output/qc_panel_239.json").read_text())
gp = qc_data.get("panel_genes", qc_data.get("genes", []))
if isinstance(gp[0], dict): gp = [g["gene"] if isinstance(g, dict) else g for g in gp]
panel_genes = [g for g in gp if g in tpm_df.index]
X_panel = np.zeros((len(ch), len(panel_genes)), dtype=np.float32)
for i, g in enumerate(panel_genes):
    X_panel[:, i] = Xr[:, list(tpm_df.index).index(g)]

print(f"Data: {len(ch)} samples, {len(panel_genes)} panel genes")

# ═══════════════════════════════════════════════════════════
# 1. LITERATURE BENCHMARK
# ═══════════════════════════════════════════════════════════
print("\n─── 1. Literature Benchmark ───")

# Curated gene sets from MSigDB, PanglaoDB, muscle stem cell literature
gene_sets = {
    "MSigDB_Hallmark_Myogenesis": ["ACTA1","ACTA2","ACTN2","MYH11","MYH2","MYL1","MYL9","MYLPF","TNNI1","TNNT1","TNNT2","TPM1","TPM2","DES","MYLK","MYL6","MYL6B","MYOM1","TTN","NEB","CSRP3","LDB3","MYOZ1","TCAP","ANKRD1","CKM","CKMT2","MB","PGAM2","PYGM"],
    "MSigDB_EMT": ["VIM","FN1","CDH2","SNAI1","SNAI2","TWIST1","ZEB1","ZEB2","MMP2","MMP9","COL1A1","COL1A2","COL3A1","COL5A2","ITGA5","ITGB1","SPARC","TGFB1","WNT5A","CTNNB1"],
    "PanglaoDB_Satellite_Cell": ["PAX7","MYF5","MYOD1","MYOG","CD34","CXCR4","CALCR","VCAM1","ITGA7","CHRNA1","MEGF10","NOTCH1","NOTCH2","NOTCH3","DLL1","JAG1","HES1","HEY1","HEYL","RBPJ"],
    "PanglaoDB_Myoblast": ["MYOD1","MYF5","MYOG","MYF6","PAX3","PAX7","MEF2C","MEF2D","DES","ACTA1","CKM","MB","TNNT2","MYL1","MYH3","CDH15","ITGA7","NCAM1","MET","CXCR4"],
    "PanglaoDB_Myocyte": ["MYH1","MYH2","MYH4","MYH7","ACTA1","TNNT3","TNNC1","TNNI2","MYL1","MYL3","CKM","MB","DES","NEB","TTN","MYOM1","MYOM2","MYBPC1","MYBPC2","ATP2A1"],
    "PanglaoDB_Fibroblast": ["COL1A1","COL1A2","COL3A1","COL5A1","COL6A1","FN1","VIM","DCN","LUM","FAP","ACTA2","TAGLN","S100A4","FSP1","PDGFRA","PDGFRB","THY1","CD34","NT5E","ENG"],
    "Muscle_Stem_Cell_Quiescence": ["PAX7","CALCR","CD34","CXCR4","VCAM1","SPRY1","NOTCH1","NOTCH3","RBPJ","HES1","HEY1","FOXO3","CDKN1B","CDKN1C","RB1","RBL2","EZH2","BMI1","SP1","MYC"],
    "Muscle_Differentiation_Canonical": ["MYOD1","MYF5","MYOG","MYF6","MEF2A","MEF2C","MEF2D","SRF","TEAD1","TEAD4","YAP1","WWTR1","SMAD2","SMAD3","SMAD4","HDAC4","HDAC5","HDAC9","SIX1","EYA1"],
    "Metabolic_Stemness": ["HK2","PKM","LDHA","GAPDH","ENO1","PGK1","PFKL","SLC2A1","SLC16A3","MYC","HIF1A","MTOR","AKT1","PIK3CA","PTEN","AMPK","SIRT1","PGC1A","TFAM","NRF1"],
    "Senescence_Markers": ["CDKN1A","CDKN2A","CDKN2B","TP53","RB1","GLB1","LMNB1","HMGB1","IL6","IL8","CXCL8","CCL2","MMP3","SERPINE1","IGFBP3","IGFBP5","IGFBP7","TNF","NFKB1","RELA"],
}

# Compute overlap
benchmark_results = []
for set_name, genes in sorted(gene_sets.items()):
    overlap = sorted(set(genes) & set(panel_genes))
    jaccard = len(overlap) / len(set(genes) | set(panel_genes))
    benchmark_results.append({
        "gene_set": set_name, "size": len(genes),
        "overlap_count": len(overlap), "overlap_genes": overlap,
        "jaccard": round(jaccard, 4),
        "pct_of_panel": round(100*len(overlap)/len(panel_genes), 1),
        "pct_of_set": round(100*len(overlap)/len(genes), 1)
    })
    print(f"  {set_name}: {len(overlap)}/{len(genes)} overlap, Jaccard={jaccard:.3f} -> {overlap}")

# ═══════════════════════════════════════════════════════════
# 2. DRUG/COMPOUND PREDICTION
# ═══════════════════════════════════════════════════════════
print("\n─── 2. Drug/Compound Prediction ───")

# Known compounds affecting muscle stem cell states (literature-curated)
compounds = {
    "p38_inhibitor_SB203580": {"target": "MAPK14", "effect": "maintains_quiescence", "ref": "Bernet 2014"},
    "p38_inhibitor_SB202190": {"target": "MAPK14", "effect": "promotes_self_renewal", "ref": "Cosgrove 2014"},
    "STAT3_inhibitor_Stattic": {"target": "STAT3", "effect": "promotes_expansion", "ref": "Tierney 2014"},
    "Notch_activator_Jagged1": {"target": "NOTCH1", "effect": "maintains_stemness", "ref": "Conboy 2003"},
    "Wnt_activator_CHIR99021": {"target": "GSK3B", "effect": "promotes_expansion", "ref": "Sato 2016"},
    "TGFb_inhibitor_SB431542": {"target": "TGFBR1", "effect": "blocks_differentiation", "ref": "Hicks 2018"},
    "BMP_inhibitor_Noggin": {"target": "BMP2/4", "effect": "maintains_undifferentiated", "ref": "Ono 2011"},
    "FGF2_basic_FGF": {"target": "FGFR1", "effect": "promotes_proliferation", "ref": "Allen 1999"},
    "IGF1": {"target": "IGF1R", "effect": "promotes_hypertrophy", "ref": "Musaro 2001"},
    "HDAC_inhibitor_TrichostatinA": {"target": "HDAC1/2", "effect": "promotes_myogenic", "ref": "Iezzi 2004"},
    "ROCK_inhibitor_Y27632": {"target": "ROCK1/2", "effect": "promotes_survival", "ref": "Watanabe 2007"},
    "mTOR_inhibitor_Rapamycin": {"target": "MTOR", "effect": "promotes_quiescence", "ref": "Rodgers 2014"},
    "AMPK_activator_AICAR": {"target": "PRKAA1", "effect": "metabolic_reprogramming", "ref": "Canto 2010"},
    "SIRT1_activator_Resveratrol": {"target": "SIRT1", "effect": "mitochondrial_biogenesis", "ref": "Lagouge 2006"},
    "PPARg_agonist_Rosiglitazone": {"target": "PPARG", "effect": "adipogenic_avoid", "ref": "Lehmann 1995"},
    "DNA_methylation_5azacytidine": {"target": "DNMT1", "effect": "epigenetic_reprogramming", "ref": "Taylor 1984"},
    "LSD1_inhibitor_Tranylcypromine": {"target": "KDM1A", "effect": "promotes_myogenic", "ref": "Tosic 2018"},
    "EZH2_inhibitor_GSK126": {"target": "EZH2", "effect": "promotes_differentiation", "ref": "Juan 2011"},
    "Betaine_trimethylglycine": {"target": "BHMT", "effect": "osmolyte_myogenic", "ref": "Senesi 2016"},
    "Creatine": {"target": "CKM", "effect": "energy_metabolism", "ref": "Deldicque 2008"},
}

# Score compounds by relevance to our panel genes
compound_scores = []
for name, info in compounds.items():
    target = info["target"]
    # Check if target or related genes are in our panel
    related_genes = []
    for g in panel_genes:
        if target.upper() in g.upper() or g.upper() in target.upper():
            related_genes.append(g)
    score = len(related_genes)
    compound_scores.append({
        "compound": name, "target": target, "effect": info["effect"],
        "ref": info["ref"], "panel_gene_matches": related_genes,
        "relevance_score": score
    })

compound_scores.sort(key=lambda x: x["relevance_score"], reverse=True)
print("Top compounds for media optimization:")
for c in compound_scores[:10]:
    print(f"  {c['compound']}: {c['effect']} (target={c['target']}, matches={c['panel_gene_matches']})")

# ═══════════════════════════════════════════════════════════
# 3. MINIMAL PANEL ANALYSIS
# ═══════════════════════════════════════════════════════════
print("\n─── 3. Minimal Panel Analysis ───")

# Sequential backward feature elimination
lr = LogisticRegression(C=1.0, max_iter=2000, random_state=42)
n_features_range = range(5, len(panel_genes), 5)
minimal_results = []

for n_feat in n_features_range:
    sfs = SequentialFeatureSelector(lr, n_features_to_select=n_feat, direction='backward', cv=3, n_jobs=-1)
    sfs.fit(X_panel, y)
    selected_idx = np.where(sfs.support_)[0]
    selected_genes = [panel_genes[i] for i in selected_idx]
    scores = cross_val_score(lr, X_panel[:, selected_idx], y, cv=5)
    minimal_results.append({
        "n_genes": n_feat, "genes": selected_genes,
        "mean_accuracy": round(float(scores.mean()), 4),
        "std_accuracy": round(float(scores.std()), 4)
    })
    print(f"  {n_feat} genes: {scores.mean():.3f} +/- {scores.std():.3f}")

# Find knee point (where adding more genes gives diminishing returns)
accs = [r["mean_accuracy"] for r in minimal_results]
deltas = [accs[i] - accs[i-1] for i in range(1, len(accs))]
knee_idx = next((i for i, d in enumerate(deltas) if d < 0.005), len(deltas) - 1)
knee_n = list(n_features_range)[knee_idx + 1]
print(f"\nDiminishing returns after {knee_n} genes (delta < 0.005)")

# ═══════════════════════════════════════════════════════════
# 4. qPCR NOISE SIMULATION
# ═══════════════════════════════════════════════════════════
print("\n─── 4. qPCR Noise Simulation ───")

rng = np.random.RandomState(42)
noise_levels = [0.05, 0.10, 0.15, 0.20, 0.25]  # CV levels
noise_results = []

for cv in noise_levels:
    sim_scores = []
    for _ in range(200):
        noise = rng.normal(0, cv * X_panel.std(axis=0), X_panel.shape)
        X_noisy = X_panel + noise
        X_noisy = np.maximum(X_noisy, 0)
        score = cross_val_score(lr, X_noisy, y, cv=5).mean()
        sim_scores.append(score)
    sim_scores = np.array(sim_scores)
    noise_results.append({
        "cv_level": cv, "mean_accuracy": round(float(sim_scores.mean()), 4),
        "std_accuracy": round(float(sim_scores.std()), 4),
        "min_accuracy": round(float(sim_scores.min()), 4),
        "accuracy_loss": round(0.967 - float(sim_scores.mean()), 4)
    })
    print(f"  CV={cv:.2f}: accuracy={sim_scores.mean():.3f}+/-{sim_scores.std():.3f}, loss={0.967-sim_scores.mean():.4f}")

# ═══════════════════════════════════════════════════════════
# SAVE
# ═══════════════════════════════════════════════════════════
results = {
    "literature_benchmark": benchmark_results,
    "drug_predictions": compound_scores,
    "minimal_panel": {"results": minimal_results, "knee_point_n_genes": knee_n},
    "qpcr_noise_simulation": noise_results
}
json.dump(results, open(OUT / "literature_drug_panel_noise.json", "w"), indent=2)
print("\nSaved to literature_drug_panel_noise.json")
print("DONE")
