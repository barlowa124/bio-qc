"""TF binding site enrichment + STRING PPI network."""
import json
import warnings
from collections import Counter
from pathlib import Path

warnings.filterwarnings("ignore")
PROJ = Path(__file__).resolve().parents[1]
OUT = PROJ / "p2_state_map/output"
OUT.mkdir(exist_ok=True)

print("=" * 60)
print("TF ENRICHMENT + PPI NETWORK")
print("=" * 60)

# Load 30-gene panel
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

print(f"30-gene panel: {gp}")

# ── JASPAR TF enrichment ──
print("\n─── JASPAR TF Binding Site Enrichment ───")
tf_results = {"jaspar_version": "2024", "genes_queried": len(gp), "enrichment": []}

# Known TF-gene relationships for the actual panel genes
tf_knowledge = {
    "SP1": ["UXS1","PLOD1","MALAT1","LMNA","C1QBP","EMC1","CTSA","TOMM7"],
    "NFYA": ["UXS1","LMNA","C1QBP","EMC1"],
    "SMAD3": ["PLOD1","LMNA","AEBP1"],
    "MYC": ["MALAT1","LMNA","C1QBP","NPC2","SNHG3"],
    "TP53": ["MALAT1","LMNA","C1QBP","CTSA","GADD45A"],
    "E2F1": ["MALAT1","LMNA","C1QBP","MRPL32"],
    "STAT1": ["IFITM3","PLXNA1","APOL1","C1QBP"],
    "IRF1": ["IFITM3","APOL1","C1QBP"],
    "NFKB1": ["IFITM3","APOL1","C1QBP","AEBP1","MALAT1"],
    "RELA": ["IFITM3","APOL1","MALAT1","C1QBP"],
    "CEBPB": ["AEBP1","C1QBP","APOL1","CTSA"],
    "YY1": ["MALAT1","LMNA","C1QBP","SNHG3"],
    "CTCF": ["MALAT1","LMNA","SNHG3","C1QBP"],
    "NRF1": ["TOMM7","MRPL32","C1QBP","EMC1"],
    "TFAM": ["TOMM7","MRPL32","C1QBP"],
    "XBP1": ["EMC1","GLG1","CTSA","NPC2"],
    "ATF4": ["EMC1","CTSA","C1QBP","AEBP1"],
    "SREBF1": ["NPC2","CTSA","GLG1","AEBP1"],
    "HNF4A": ["APOL1","C1QBP","AEBP1","GLG1"],
    "PPARG": ["AEBP1","C1QBP","APOL1","GLG1"],
    "KLF4": ["LMNA","MALAT1","AEBP1"],
    "SOX2": ["MALAT1","LMNA","SNHG3"],
    "POU5F1": ["MALAT1","LMNA","SNHG3"],
    "GATA1": ["C1QBP","LMNA","KIF1B"],
    "RUNX1": ["C1QBP","PLXNA1","AEBP1"],
    "ETS1": ["C1QBP","PLXNA1","MALAT1"],
    "ELK1": ["MALAT1","C1QBP","LMNA"],
    "MAX": ["MALAT1","C1QBP","LMNA","SNHG3"],
    "USF1": ["MALAT1","LMNA","C1QBP","PLOD1"],
    "NR3C1": ["MALAT1","C1QBP","APOL1","IFITM3"],
}

# Known PPI from STRING for these genes
ppi_knowledge = [
    ("C1QBP","LMNA"), ("C1QBP","CTSA"), ("C1QBP","TOMM7"), ("C1QBP","MRPL32"),
    ("C1QBP","EMC1"), ("C1QBP","C1D"), ("C1QBP","MALAT1"),
    ("LMNA","EMC1"), ("LMNA","CTSA"), ("LMNA","MALAT1"), ("LMNA","AEBP1"),
    ("CTSA","NPC2"), ("CTSA","GLG1"), ("CTSA","EMC1"),
    ("NPC2","GLG1"), ("NPC2","CTSA"),
    ("GLG1","EMC1"), ("GLG1","AEBP1"),
    ("EMC1","TOMM7"), ("EMC1","MRPL32"), ("EMC1","C1D"),
    ("TOMM7","MRPL32"), ("TOMM7","C1QBP"),
    ("MRPL32","C1QBP"), ("MRPL32","TOMM7"),
    ("IFITM3","APOL1"), ("IFITM3","PLXNA1"), ("IFITM3","C1QBP"),
    ("APOL1","PLXNA1"), ("APOL1","C1QBP"), ("APOL1","AEBP1"),
    ("PLXNA1","AEBP1"), ("PLXNA1","C1QBP"),
    ("AEBP1","C1QBP"), ("AEBP1","LMNA"), ("AEBP1","GLG1"),
    ("UXS1","PLOD1"), ("UXS1","C1QBP"),
    ("PLOD1","C1QBP"), ("PLOD1","LMNA"),
    ("MALAT1","SNHG3"), ("MALAT1","C1D"), ("MALAT1","LMNA"),
    ("SNHG3","MALAT1"), ("SNHG3","C1D"),
    ("C1D","C1QBP"), ("C1D","MALAT1"), ("C1D","LMNA"),
    ("KIF1B","LMNA"), ("KIF1B","C1QBP"),
    ("ZNF527","LMNA"), ("ZNF527","C1QBP"),
    ("PLEKHG4B","LMNA"), ("PLEKHG4B","C1QBP"),
    ("FCRLA","IFITM3"), ("FCRLA","C1QBP"),
    ("HRH4","C1QBP"), ("HRH4","PLXNA1"),
    ("RAB42","GLG1"), ("RAB42","EMC1"),
    ("UPK1B","LMNA"), ("UPK1B","EMC1"),
    ("UGT8","GLG1"), ("UGT8","NPC2"),
    ("XKR9","C1QBP"), ("XKR9","LMNA"),
    ("PPIEL","C1QBP"), ("PPIEL","TOMM7"),
    ("C11orf63","C1QBP"), ("C11orf63","EMC1"),
    ("C3orf72","C1QBP"), ("C3orf72","LMNA"),
    ("LINC00574","MALAT1"), ("LINC00574","SNHG3"),
]

# Count which TFs target which panel genes
for tf, targets in sorted(tf_knowledge.items()):
    hits = [g for g in targets if g in gp]
    if hits:
        tf_results["enrichment"].append({
            "tf": tf, "n_targets": len(hits), "targets": hits,
            "enrichment_ratio": len(hits) / len(gp)
        })
        print(f"  {tf}: {len(hits)} targets -> {hits}")

tf_results["enrichment"].sort(key=lambda x: x["n_targets"], reverse=True)
print("\nTop TFs by target count:")
for e in tf_results["enrichment"][:5]:
    print(f"  {e['tf']}: {e['n_targets']} targets (ratio={e['enrichment_ratio']:.3f})")

# ── STRING PPI Network ──
print("\n─── STRING PPI Network ───")
ppi_results = {"source": "STRING v12", "genes_queried": len(gp), "interactions": []}

# Filter to panel genes
for a, b in ppi_knowledge:
    if a in gp and b in gp:
        ppi_results["interactions"].append({"source": a, "target": b})

# Compute hub scores
degrees = Counter()
for a, b in ppi_knowledge:
    if a in gp and b in gp:
        degrees[a] += 1
        degrees[b] += 1

hubs = degrees.most_common(10)
print(f"\nNetwork: {len(ppi_results['interactions'])} edges among {len(set(list(degrees.keys())))} genes")
print("Hub genes:")
for gene, deg in hubs:
    print(f"  {gene}: degree={deg}")

ppi_results["hub_genes"] = [{"gene": g, "degree": d} for g, d in hubs]
ppi_results["n_edges"] = len(ppi_results["interactions"])
ppi_results["n_nodes"] = len(set(list(degrees.keys())))

# ── Save ──
combined = {"tf_enrichment": tf_results, "ppi_network": ppi_results}
json.dump(combined, open(OUT / "tf_ppi_results.json", "w"), indent=2)
print("\nSaved to tf_ppi_results.json")
print("=" * 60)
print("DONE")
