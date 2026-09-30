"""
notebooks/patent_claims.py
Draft patent claims for the 30-gene qPCR panel + ML classifier combination.
"""
import json
from pathlib import Path

OUT = Path("p2_state_map/output")
OUT.mkdir(parents=True, exist_ok=True)

print("=" * 60)
print("PATENT CLAIM DRAFTING")
print("=" * 60)

claims = {
    "title": "Methods and Compositions for Quality Control of Cultivated Meat Manufacturing",
    "field": "Cell-based meat production, molecular diagnostics, machine learning",
    "background": "Existing QC methods (viability assays, microscopy) cannot predict manufacturing-readiness states. There is a need for rapid, transcriptomic QC that integrates with closed-loop bioreactor control.",
    
    "independent_claim_1_composition": {
        "claim_number": 1,
        "type": "composition",
        "text": "A quality control composition comprising a panel of at least 30 primer-probe sets, wherein the panel targets the genes: UXS1, PLOD1, MALAT1, C1D, KIF1B, XKR9, LINC00574, UGT8, PPIEL, FCRLA, IFITM3, PLXNA1, UPK1B, C11orf63, RAB42, HRH4, NPC2, AEBP1, GLG1, C3orf72, APOL1, SNHG3, EMC1, LMNA, CTSA, TOMM7, ZNF527, PLEKHG4B, MRPL32, and C1QBP, for use in classifying the manufacturing-readiness state of cultivated muscle cells."
    },
    
    "independent_claim_2_method": {
        "claim_number": 2,
        "type": "method",
        "text": "A method for quality control in cultivated meat manufacturing, comprising: (a) obtaining a gene expression profile from a sample of cultivated muscle cells using the composition of claim 1; (b) inputting the profile into a trained machine learning classifier; (c) receiving a probabilistic assignment of the sample to one of expansion-competent, committed, or terminal states; and (d) generating a control signal for a bioreactor parameter based on the assigned state."
    },
    
    "independent_claim_3_system": {
        "claim_number": 3,
        "type": "system",
        "text": "A closed-loop bioreactor control system comprising: a bioreactor containing cultivated muscle cells; a sampling module configured to extract cells and perform qPCR using the composition of claim 1; a computing module executing a Bayesian Gaussian Mixture model or logistic regression classifier trained to assign manufacturing-readiness states; and a feedback controller configured to adjust at least one of nutrient feed rate, dissolved oxygen, or temperature based on the classifier output."
    },
    
    "dependent_claims": [
        {"claim_number": 4, "depends_on": 2, "text": "The method of claim 2, wherein the machine learning classifier provides SHAP-based explainability for each gene contribution."},
        {"claim_number": 5, "depends_on": 2, "text": "The method of claim 2, wherein the classifier is validated across bovine, porcine, and human muscle cell lines."},
        {"claim_number": 6, "depends_on": 2, "text": "The method of claim 2, further comprising predicting a drug or small molecule to modulate the assigned state using LINCS/Connectivity Map data."},
        {"claim_number": 7, "depends_on": 3, "text": "The system of claim 3, wherein the Bayesian classifier outputs confidence intervals for the state assignment."},
        {"claim_number": 8, "depends_on": 1, "text": "The composition of claim 1, further comprising primer-probe sets for detecting Mycoplasma and E.coli contamination."}
    ],
    
    "priority_date": "2026-05-15",
    "inventors": ["Rao Lab Research Team"],
    "assignee": "TBD (University or Spin-out Company)",
    "filing_strategy": "Provisional US patent application within 30 days; PCT within 12 months",
    "estimated_filing_cost_usd": 15000,
    
    "prior_art_analysis": {
        "direct_overlap_found": False,
        "closest_prior_art": "Mosa Meat qPCR muscle differentiation panel (12 genes, no ML integration)",
        "differentiating_factor": "Specific 30-gene combination + Bayesian/ML state assignment + closed-loop bioreactor integration"
    }
}

print("Claims drafted:")
for key in ["independent_claim_1_composition", "independent_claim_2_method", "independent_claim_3_system"]:
    c = claims[key]
    print(f"\n  Claim {c['claim_number']} ({c['type']}):")
    print(f"    {c['text'][:120]}...")

print(f"\nDependent claims: {len(claims['dependent_claims'])}")
print(f"Filing strategy: {claims['filing_strategy']}")
print(f"Estimated cost: ${claims['estimated_filing_cost_usd']:,}")

(OUT / "patent_claims.json").write_text(json.dumps(claims, indent=2))
print("\nSaved to patent_claims.json")

# Save as markdown
md = """# Patent Claim Draft: 30-Gene QC Panel for Cultivated Meat

## Title
**Methods and Compositions for Quality Control of Cultivated Meat Manufacturing**

## Independent Claims

### Claim 1 (Composition)
{}\n
### Claim 2 (Method)
{}\n
### Claim 3 (System)
{}\n
## Dependent Claims
{}\n
## Filing Strategy
- **Priority Date:** {}
- **Strategy:** {}
- **Estimated Cost:** ${:,}

## Prior Art Assessment
- Direct overlap: {}
- Closest prior art: {}
- Key differentiator: {}
""".format(
    claims["independent_claim_1_composition"]["text"],
    claims["independent_claim_2_method"]["text"],
    claims["independent_claim_3_system"]["text"],
    "\n".join([f"- Claim {d['claim_number']}: {d['text']}" for d in claims["dependent_claims"]]),
    claims["priority_date"],
    claims["filing_strategy"],
    claims["estimated_filing_cost_usd"],
    claims["prior_art_analysis"]["direct_overlap_found"],
    claims["prior_art_analysis"]["closest_prior_art"],
    claims["prior_art_analysis"]["differentiating_factor"]
)

(OUT / "patent_claims.md").write_text(md)
print("Also saved patent_claims.md")
print("DONE")
