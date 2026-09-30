"""Deepen P1: Extract protocol factors from docx text + P4: quantify media components.
Plus cross-project integration, METAFlux mapping, and validation design."""
import json
import re
import warnings
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from docx import Document

warnings.filterwarnings("ignore")
PROJ = Path(__file__).resolve().parents[1]
OUT1 = PROJ / "p1_protocol_meta_analysis/output"
OUT4 = PROJ / "p4_media_formulation/output"
OUTX = PROJ / "p2_state_map/output"
OUT1.mkdir(exist_ok=True); OUT4.mkdir(exist_ok=True)

print("=" * 60)
print("DEEPENED P1 + P4 + INTEGRATION + MAPPING + VALIDATION")
print("=" * 60)

# ═══════════════════════════════════════════════════════════
# P1 DEEPENED: Extract protocol factors
# ═══════════════════════════════════════════════════════════
print("\n─── P1: Deep Protocol Factor Extraction ───")

protocol_dir = Path(r"D:\Rao_Lab_Google_Drive_Downloads\Stem Cell Lab\RJ Taylor\Protocols")
main_protocols = r"D:\Rao_Lab_Google_Drive_Downloads\Stem Cell Lab\Protocols.docx"

all_text = ""
protocol_sections = {}

# Parse main Protocols.docx
try:
    doc = Document(main_protocols)
    for p in doc.paragraphs:
        if p.text.strip():
            all_text += p.text + "\n"
            # Detect section headers
            if p.style.name.startswith("Heading") or p.text.strip().isupper() or len(p.text.strip()) < 80 and p.text.strip().endswith(":"):
                protocol_sections[p.text.strip()] = []
    print(f"   Main Protocols.docx: {len(doc.paragraphs)} paragraphs, {len(protocol_sections)} sections")
except Exception as e:
    print(f"   Main Protocols.docx error: {e}")

# Parse RJ Taylor protocols
rj_texts = {}
for f in protocol_dir.glob("*.docx"):
    try:
        doc = Document(str(f))
        text = "\n".join([p.text for p in doc.paragraphs if p.text.strip()])
        rj_texts[f.stem] = text
        all_text += text + "\n"
        print(f"   {f.name}: {len(text)} chars")
    except Exception as e:
        print(f"   {f.name}: {e}")

# Extract protocol factors
factors = {
    "media_components": [],
    "timing": [],
    "equipment": [],
    "cell_types": [],
    "concentrations": [],
    "steps": [],
}

# Media components
media_keywords = ["DMEM", "RPMI", "FBS", "BSA", "serum", "antibiotic", "penicillin", "streptomycin", "glutamine", "glucose", "sodium pyruvate", "NEAA", "bFGF", "EGF", "TGF", "insulin", "transferrin", "selenium", "hydrocortisone", "ascorbic acid", "2-mercaptoethanol", "DMSO", "collagenase", "trypsin", "accutase", "dispase", "Y-27632", "ROCK inhibitor", "CHIR99021", "PD0325901", "LIF", "BMP4", "Activin A", "Wnt", "FGF2", "FGF4", "VEGF", "PDGF", "IGF"]
for kw in media_keywords:
    count = len(re.findall(rf"\b{re.escape(kw)}\b", all_text, re.IGNORECASE))
    if count > 0:
        factors["media_components"].append({"component": kw, "mentions": count})

factors["media_components"].sort(key=lambda x: x["mentions"], reverse=True)

# Timing
time_patterns = [r"(\d+)[\s-]*hour", r"(\d+)[\s-]*minute", r"(\d+)[\s-]*day", r"(\d+)[\s-]*week"]
for pat, unit in zip(time_patterns, ["hours", "minutes", "days", "weeks"]):
    matches = re.findall(pat, all_text, re.IGNORECASE)
    if matches:
        factors["timing"].append({"unit": unit, "values": [int(m) for m in matches[:20]], "count": len(matches)})

# Concentrations
conc_pattern = r"(\d+\.?\d*)\s*(μg/ml|ng/ml|mg/ml|μM|nM|mM|%|U/ml|IU/ml)"
conc_matches = re.findall(conc_pattern, all_text, re.IGNORECASE)
conc_counter = Counter([f"{v} {u}" for v, u in conc_matches])
factors["concentrations"] = [{"value": k, "count": v} for k, v in conc_counter.most_common(20)]

# Steps
step_patterns = [r"(\d+)\s*\.\s+([A-Z][^\n]+)", r"Step\s+(\d+)[:\s]+([^\n]+)", r"([A-Z][^\n]+(?:incubate|wash|centrifuge|plate|seed|feed|change|passage|freeze|thaw|fix|stain|extract|purify|measure)[^\n]*)"] 
for pat in step_patterns:
    matches = re.findall(pat, all_text, re.IGNORECASE)
    for m in matches[:30]:
        if isinstance(m, tuple):
            factors["steps"].append(" ".join(m)[:120])
        else:
            factors["steps"].append(m[:120])

print(f"   Media components found: {len(factors['media_components'])}")
print(f"   Timing patterns: {len(factors['timing'])}")
print(f"   Concentration values: {len(factors['concentrations'])}")
print(f"   Steps extracted: {len(factors['steps'])}")

# ═══════════════════════════════════════════════════════════
# P4 DEEPENED: Quantify media components
# ═══════════════════════════════════════════════════════════
print("\n─── P4: Deep Media Component Quantification ───")

media_analysis = {}

# Parse Media.xlsx
try:
    media = pd.read_excel(r"D:\Rao_Lab_Google_Drive_Downloads\Stem Cell Lab\Media.xlsx", header=None)
    # Try to find structured data
    media_analysis["media_shape"] = list(media.shape)
    # Extract any numeric columns
    numeric_cols = media.select_dtypes(include=[np.number]).columns.tolist()
    media_analysis["numeric_columns"] = len(numeric_cols)
    # Look for component names in first column
    first_col = media.iloc[:, 0].dropna().astype(str).tolist()
    component_like = [v for v in first_col if len(v) > 3 and not v.startswith("Unnamed")][:30]
    media_analysis["potential_components"] = component_like
    print(f"   Media.xlsx: {media.shape}, {len(numeric_cols)} numeric cols, {len(component_like)} potential components")
except Exception as e:
    print(f"   Media.xlsx error: {e}")

# Parse Ordering.xlsx
try:
    ordering = pd.read_excel(r"D:\Rao_Lab_Google_Drive_Downloads\Stem Cell Lab\Ordering.xlsx")
    media_analysis["ordering_shape"] = list(ordering.shape)
    # Look for cost data
    cost_cols = [c for c in ordering.columns if any(k in str(c).lower() for k in ["cost", "price", "$", "amount"])]
    media_analysis["cost_columns"] = cost_cols
    print(f"   Ordering.xlsx: {ordering.shape}, cost cols: {cost_cols}")
except Exception as e:
    print(f"   Ordering.xlsx error: {e}")

# Parse CEF Media Experiments
try:
    cef = pd.read_excel(r"D:\Rao_Lab_Google_Drive_Downloads\Stem Cell Lab\Angus Barlow\Copy of CEF Media Experiments.xlsx")
    media_analysis["cef_shape"] = list(cef.shape)
    media_analysis["cef_columns"] = cef.columns.tolist()
    print(f"   CEF Media: {cef.shape}, cols: {cef.columns.tolist()}")
except Exception as e:
    print(f"   CEF Media error: {e}")

# ═══════════════════════════════════════════════════════════
# CROSS-PROJECT INTEGRATION
# ═══════════════════════════════════════════════════════════
print("\n─── Cross-Project Integration ───")

# Load state map data
DATA = Path(r"C:\Users\asdf\CascadeProjects\rao_lab_ml\nn_pipeline\combined_data\metaflux_scfea_pseudobulk")
tpm = pd.read_csv(DATA / "scfea_pseudobulk_tpm.csv", index_col=0)
flux = pd.read_csv(DATA / "metaflux_combined_flux_matrix_all.csv", index_col=0)
info = pd.read_csv(DATA / "sample_info.csv")
common = sorted(set(tpm.columns) & set(flux.columns))

X_rna = np.log1p(tpm[common].values.T).astype(np.float32)
X_flux = flux[common].values.T.astype(np.float32)
genes = tpm.index.values

gvar = X_rna.var(axis=0)
keep = gvar > np.percentile(gvar, 25)
X_rna_f = X_rna[:, keep]; kept = genes[keep]

from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

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

# Integration: per-dataset readiness distribution
sample_datasets = []
for s in common:
    row = info[info["sample_id"] == s]
    ds = row["dataset"].values[0] if len(row) > 0 else "unknown"
    sample_datasets.append(ds)

datasets = sorted(set(sample_datasets))
integration = {}
for ds in datasets:
    idx = [i for i, s in enumerate(common) if sample_datasets[i] == ds]
    if len(idx) >= 3:
        dist = Counter(readiness[idx])
        integration[ds] = {"n": len(idx), "expansion_pct": dist.get("expansion_competent", 0)/len(idx), "committed_pct": dist.get("committed", 0)/len(idx), "terminal_pct": dist.get("terminal", 0)/len(idx)}
        print(f"   {ds}: n={len(idx)}, expansion={integration[ds]['expansion_pct']:.0%}, committed={integration[ds]['committed_pct']:.0%}, terminal={integration[ds]['terminal_pct']:.0%}")

# ═══════════════════════════════════════════════════════════
# METAFLUX REACTION-TO-PATHWAY MAPPING
# ═══════════════════════════════════════════════════════════
print("\n─── METAFlux Reaction-to-Pathway Mapping ───")

# Map HMR reactions to 10 PANN pathways based on reaction names
pathway_keywords = {
    "OXPHOS": ["oxphos", "oxidative phosphorylation", "complex I", "complex II", "complex III", "complex IV", "ATP synthase", "cytochrome", "ubiquinone", "NADH dehydrogenase", "succinate dehydrogenase", "coenzyme Q", "electron transport", "HMR_6916", "HMR_6917", "HMR_6918", "HMR_6919", "HMR_6920", "HMR_6921", "HMR_6922", "HMR_6923", "HMR_6924", "HMR_6925"],
    "Glycolysis": ["glycolysis", "glucokinase", "hexokinase", "phosphoglucose", "phosphofructokinase", "aldolase", "triose", "glyceraldehyde", "phosphoglycerate", "enolase", "pyruvate kinase", "lactate dehydrogenase", "glucose-6-phosphate", "fructose-6-phosphate", "fructose-1,6-bisphosphate", "DHAP", "G3P", "1,3-BPG", "3-PG", "2-PG", "PEP", "pyruvate", "lactate"],
    "TCA_Cycle": ["TCA", "citrate", "isocitrate", "alpha-ketoglutarate", "succinyl-CoA", "succinate", "fumarate", "malate", "oxaloacetate", "aconitase", "citrate synthase", "isocitrate dehydrogenase", "alpha-ketoglutarate dehydrogenase", "succinyl-CoA synthetase", "succinate dehydrogenase", "fumarase", "malate dehydrogenase"],
    "PPP": ["pentose phosphate", "glucose-6-phosphate dehydrogenase", "6-phosphogluconate", "ribulose-5-phosphate", "ribose-5-phosphate", "transketolase", "transaldolase", "xylulose-5-phosphate", "sedoheptulose-7-phosphate", "erythrose-4-phosphate", "G6PD", "6PGD", "NADPH"],
    "Fatty_Acid_Oxidation": ["fatty acid oxidation", "beta-oxidation", "carnitine", "acyl-CoA", "acetyl-CoA", "palmitoyl", "stearoyl", "oleoyl", "linoleoyl", "fatty acyl", "acylcarnitine", "CPT1", "CPT2", "acyl-CoA dehydrogenase", "enoyl-CoA", "hydroxyacyl-CoA", "ketoacyl-CoA", "thiolase"],
    "Amino_Acid_Metabolism": ["amino acid", "glutamate", "glutamine", "aspartate", "asparagine", "alanine", "serine", "glycine", "cysteine", "methionine", "valine", "leucine", "isoleucine", "phenylalanine", "tyrosine", "tryptophan", "histidine", "lysine", "arginine", "proline", "threonine", "transaminase", "aminotransferase", "deaminase", "decarboxylase"],
    "Nucleotide_Metabolism": ["nucleotide", "purine", "pyrimidine", "adenine", "guanine", "cytosine", "thymine", "uracil", "adenosine", "guanosine", "cytidine", "thymidine", "uridine", "AMP", "GMP", "CMP", "TMP", "UMP", "ATP", "GTP", "CTP", "TTP", "UTP", "IMP", "XMP", "ribonucleotide", "deoxyribonucleotide"],
    "Lipid_Metabolism": ["lipid", "cholesterol", "sterol", "phospholipid", "sphingolipid", "ceramide", "triglyceride", "triacylglycerol", "diacylglycerol", "monoacylglycerol", "phosphatidyl", "mevalonate", "HMG-CoA", "squalene", "lanosterol", "desmosterol", "lipoprotein", "chylomicron"],
    "ROS_Detox": ["ROS", "reactive oxygen", "superoxide", "hydrogen peroxide", "catalase", "peroxidase", "glutathione", "thioredoxin", "peroxiredoxin", "superoxide dismutase", "SOD", "GPX", "GSH", "GSSG", "NADPH oxidase", "oxidative stress", "antioxidant"],
    "One_Carbon_Metabolism": ["one-carbon", "folate", "tetrahydrofolate", "methionine", "SAM", "SAH", "homocysteine", "methyltransferase", "methylation", "formyl", "methylene", "methenyl", "MTHFR", "MTR", "BHMT", "SHMT", "serine hydroxymethyltransferase"],
}

reaction_names = flux.index.tolist()
mapping = {pw: [] for pw in pathway_keywords}
unmapped = []

for rxn in reaction_names:
    rxn_lower = rxn.lower()
    mapped = False
    for pw, keywords in pathway_keywords.items():
        if any(kw.lower() in rxn_lower for kw in keywords):
            mapping[pw].append(rxn)
            mapped = True
            break
    if not mapped:
        unmapped.append(rxn)

for pw, rxns in mapping.items():
    print(f"   {pw}: {len(rxns)} reactions")

print(f"   Unmapped: {len(unmapped)}/{len(reaction_names)} reactions")

# ═══════════════════════════════════════════════════════════
# PROSPECTIVE VALIDATION DESIGN
# ═══════════════════════════════════════════════════════════
print("\n─── Prospective Validation Design ───")

validation_design = {
    "objective": "Prospectively validate the 30-gene QC panel for batch triage in Rao Lab stem cell cultures",
    "sample_size": "Minimum 30 batches (10 expansion-competent, 10 committed, 10 terminal expected)",
    "design": "Longitudinal sampling at 3 timepoints per batch: day 0 (seeding), day 3 (mid-expansion), day 7 (pre-differentiation)",
    "measurements": [
        "30-gene qPCR panel (primary)",
        "Morphology scoring (current standard, for comparison)",
        "Viability (trypan blue)",
        "Proliferation rate (doubling time)",
        "Differentiation efficiency (post-hoc, for ground truth)",
    ],
    "endpoints": {
        "primary": "Concordance between qPCR panel prediction and actual batch outcome (pass/fail)",
        "secondary": "Time saved vs retrospective assessment, cost per batch, operator agreement",
    },
    "analysis_plan": [
        "Compare qPCR panel accuracy vs morphology scoring for batch triage",
        "ROC analysis with differentiation efficiency as ground truth",
        "Inter-operator variability assessment (3 operators, blinded)",
        "Cost-benefit analysis: qPCR ($50-100) vs failed batch cost ($X,000)",
    ],
    "timeline": "8-12 weeks for 30 batches through full expansion-differentiation cycle",
    "required_resources": [
        "qPCR machine + reagents for 30-gene panel",
        "3 trained operators for inter-operator study",
        "Cell culture consumables for 30 batches",
        "RNA extraction kit + cDNA synthesis reagents",
    ],
}

# ═══════════════════════════════════════════════════════════
# SAVE ALL
# ═══════════════════════════════════════════════════════════
print("\n─── Saving All Results ───")

json.dump({"protocol_factors": factors, "rj_texts": {k: v[:500] for k, v in rj_texts.items()}}, open(OUT1 / "protocol_factors_deep.json", "w"), indent=2, default=str)
json.dump(media_analysis, open(OUT4 / "media_components_deep.json", "w"), indent=2, default=str)
json.dump({"per_dataset_integration": integration, "pathway_mapping_summary": {pw: len(rxns) for pw, rxns in mapping.items()}, "unmapped_count": len(unmapped), "mapping": {pw: rxns[:20] for pw, rxns in mapping.items()}}, open(OUTX / "integration_and_mapping.json", "w"), indent=2)
json.dump(validation_design, open(OUTX / "validation_design.json", "w"), indent=2)

print(f"\nOutputs: {OUT1}, {OUT4}, {OUTX}")
print("=" * 60)
print("ALL DEEPENED ANALYSES COMPLETE")
print("=" * 60)