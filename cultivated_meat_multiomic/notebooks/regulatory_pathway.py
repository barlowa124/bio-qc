"""Regulatory pathway mapping for FDA/EMA cultivated meat QC panel approval."""
import json
from datetime import datetime
from pathlib import Path

PROJ = Path(__file__).resolve().parents[1]
OUT = PROJ / "p2_state_map/output"
OUT.mkdir(parents=True, exist_ok=True)

print("=" * 60)
print("REGULATORY PATHWAY MAPPING")
print("=" * 60)

# ── FDA Framework (US) ──
print("\n─── FDA (US) Framework ───")
fda = {
    "jurisdiction": "FDA Center for Food Safety and Applied Nutrition (CFSAN)",
    "framework": "Novel Food / Food Additive / GRAS",
    "current_status": "FDA issued 'no questions' letters to UPSIDE Foods (Nov 2022) and GOOD Meat (Jun 2023)",
    "qc_panel_classification": "Analytical method / quality control tool (not a food ingredient)",
    "approval_pathway": [
        "1. Pre-submission meeting with FDA Office of Food Additive Safety (optional)",
        "2. Method validation package (accuracy, precision, LOD/LOQ, specificity)",
        "3. GRAS self-determination or GRAS notice (if panel is food-contact)",
        "4. FSMA Preventive Controls: QC panel as verification tool in hazard analysis",
        "5. USDA-FSIS joint oversight for labeling if product contains cultivated cells",
    ],
    "estimated_timeline_months": 12,
    "estimated_cost_usd": 150000,
    "key_documents": [
        "FDA Guidance for Industry: Assessing the Safety of Food Ingredients",
        "FDA FSMA Final Rule for Preventive Controls for Human Food",
        "FDA GRAS Notice Inventory (for precedent)",
    ],
    "relevant_guidance_quotes": [
        "Analytical methods used for quality control must be validated for their intended use (FDA, ICH Q2(R1))",
        "Molecular methods (PCR, qPCR) are acceptable for food safety and quality verification",
    ],
}
for k, v in fda.items():
    print(f"  {k}: {v}")

# ── EFSA/EU Framework ──
print("\n─── EFSA/EU Framework ───")
efsa = {
    "jurisdiction": "European Food Safety Authority (EFSA) + EU Commission",
    "framework": "Novel Food Regulation (EU) 2015/2283",
    "current_status": "EFSA received first cultivated meat dossier (2023); no approvals yet",
    "qc_panel_classification": "Analytical method supporting Novel Food safety dossier",
    "approval_pathway": [
        "1. Pre-submission advice from EFSA (optional, 3 months)",
        "2. Novel Food application by EU member state (12-18 months)",
        "3. EFSA scientific opinion on safety and nutritional assessment",
        "4. EU Commission implementing act authorizing the Novel Food",
        "5. QC panel can be cited as part of the manufacturing process control (MPC) section",
    ],
    "estimated_timeline_months": 24,
    "estimated_cost_usd": 250000,
    "key_documents": [
        "EFSA Guidance on the preparation and presentation of an application for authorisation of a novel food",
        "EFSA Scientific Opinion on safety of cultured meat (expected 2025)",
        "EU Regulation 2015/2283 on novel foods",
    ],
    "relevant_guidance_quotes": [
        "The manufacturing process must be sufficiently described, including quality control measures (EFSA NDA Panel)",
        "Molecular characterization of the cell line and product consistency are required",
    ],
}
for k, v in efsa.items():
    print(f"  {k}: {v}")

# ── Singapore (SFA) ──
print("\n─── Singapore (SFA) Framework ───")
sfa = {
    "jurisdiction": "Singapore Food Agency (SFA)",
    "framework": "Novel Food Safety Assessment",
    "current_status": "First approval globally (Eat Just, Dec 2020); most permissive regime",
    "qc_panel_classification": "Manufacturing process control tool",
    "approval_pathway": [
        "1. Pre-submission consultation with SFA (recommended)",
        "2. Safety assessment dossier (toxicology, allergenicity, nutritional)",
        "3. Manufacturing process description with QC checkpoints",
        "4. SFA review (6-12 months)",
    ],
    "estimated_timeline_months": 9,
    "estimated_cost_usd": 80000,
    "key_documents": [
        "SFA Guidelines on Safety Assessment of Novel Foods",
        "SFA Technical Guidance for Cell-Based Meat",
    ],
}
for k, v in sfa.items():
    print(f"  {k}: {v}")

# ── Israel (IMOH) ──
print("\n─── Israel (IMOH) Framework ───")
imoh = {
    "jurisdiction": "Israel Ministry of Health (IMOH)",
    "framework": "Novel Food / New Food Ingredient",
    "current_status": "Aleph Farms approved (Jan 2024); fast-track for cultivated meat",
    "approval_pathway": [
        "1. Pre-submission meeting with IMOH Food Division",
        "2. Safety dossier (toxicology, allergenicity, stability)",
        "3. Manufacturing process with in-process controls",
        "4. IMOH review (9-15 months)",
    ],
    "estimated_timeline_months": 12,
    "estimated_cost_usd": 100000,
}
for k, v in imoh.items():
    print(f"  {k}: {v}")

# ── Regulatory milestones for the QC panel ──
print("\n─── Regulatory Milestones for 30-Gene QC Panel ───")
milestones = [
    {"phase": "Pre-clinical", "task": "Complete method validation (precision, accuracy, linearity)", "timeline": "Month 1-3", "cost": 20000},
    {"phase": "Pre-clinical", "task": "Analytical reference standard establishment", "timeline": "Month 2-4", "cost": 15000},
    {"phase": "Pre-submission", "task": "Draft FDA pre-submission meeting request", "timeline": "Month 3-4", "cost": 5000},
    {"phase": "Pre-submission", "task": "FDA pre-submission meeting", "timeline": "Month 5", "cost": 10000},
    {"phase": "Submission", "task": "Compile GRAS determination or method validation package", "timeline": "Month 6-9", "cost": 40000},
    {"phase": "Submission", "task": "Submit to FDA (US) and SFA (Singapore parallel)", "timeline": "Month 10", "cost": 15000},
    {"phase": "Review", "task": "Respond to agency questions / amendments", "timeline": "Month 11-14", "cost": 25000},
    {"phase": "Approval", "task": "FDA GRAS self-affirmation or no-objection letter", "timeline": "Month 15", "cost": 10000},
    {"phase": "Post-market", "task": "Annual surveillance / batch release records", "timeline": "Ongoing", "cost": 10000},
]
for m in milestones:
    print(f"  {m['phase']:12} | {m['timeline']:8} | ${m['cost']:>6,} | {m['task']}")

total_cost = sum(m["cost"] for m in milestones)
print(f"\n  Total regulatory cost: ${total_cost:,}")

# ── Risk register ──
print("\n─── Regulatory Risk Register ───")
risks = [
    {"risk": "FDA reclassifies qPCR panel as medical device (IVD)", "likelihood": "low", "impact": "high", "mitigation": "Frame as food QC method, not diagnostic; cite food industry qPCR precedents"},
    {"risk": "EFSA requires animal testing for Novel Food safety", "likelihood": "medium", "impact": "high", "mitigation": "Emphasize in-vitro safety data; leverage read-across from existing muscle cell cultures"},
    {"risk": "USDA-FSIS claims jurisdiction over QC method labeling", "likelihood": "low", "impact": "medium", "mitigation": "Coordinate FDA-USDA jurisdiction split early; reference Memorandum of Understanding"},
    {"risk": "Patent thickets from competitors (Mosa, Aleph)", "likelihood": "medium", "impact": "medium", "mitigation": "File provisional patent; design around existing claims by using specific gene combination"},
    {"risk": "International harmonization delays (CODEX, ISO)", "likelihood": "medium", "impact": "low", "mitigation": "Start with SFA/FDA fast-track; use mutual recognition for EU later"},
]
for r in risks:
    print(f"  [{r['likelihood'].upper()}/{r['impact'].upper()}] {r['risk']}")
    print(f"    → Mitigation: {r['mitigation']}")

# ── Save ──
results = {
    "date": datetime.now().isoformat(),
    "fda": fda,
    "efsa": efsa,
    "sfa": sfa,
    "imoh": imoh,
    "milestones": milestones,
    "total_regulatory_cost_usd": total_cost,
    "risks": risks,
    "recommendation": "Simultaneous FDA pre-submission + SFA fast-track; file provisional patent; begin method validation immediately",
}
json.dump(results, open(OUT / "regulatory_pathway.json", "w"), indent=2)
print("\nSaved to regulatory_pathway.json")
print("DONE")
