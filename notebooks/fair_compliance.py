"""FAIR compliance checklist for the cultivated meat multi-omic project."""
import json
from datetime import datetime
from pathlib import Path

PROJ = Path(__file__).resolve().parents[1]
OUT = PROJ / "p2_state_map/output"
OUT.mkdir(parents=True, exist_ok=True)

print("=" * 60)
print("FAIR COMPLIANCE CHECKLIST")
print("=" * 60)

checklist = {
    "assessment_date": datetime.now().isoformat(),
    "project": "Cultivated Meat Multi-Omic Analysis Pipeline",
    "repository": "https://github.com/barlowa124/cultivated-meat-multiomic",
    "fair_version": "1.0",
    "principles": {
        "Findable": {
            "F1_data_assigned_globally_unique_identifier": {
                "status": "partial",
                "evidence": "GitHub repo provides unique URL; DOI pending (Zenodo deposit planned)",
                "action": "Create Zenodo release to obtain DOI",
            },
            "F2_data_described_with_rich_metadata": {
                "status": "complete",
                "evidence": "README.md, FINAL_SUMMARY.md, .zenodo.json, and metadata files in output/",
                "action": "None",
            },
            "F3_metadata_clearly_explicitly_includes_identifier": {
                "status": "complete",
                "evidence": ".zenodo.json includes relatedIdentifiers linking to GitHub repo",
                "action": "None",
            },
            "F4_data_registered_indexed_searchable": {
                "status": "partial",
                "evidence": "GitHub indexed by Google; not yet in data repositories (Zenodo pending)",
                "action": "Deposit in Zenodo and Figshare; register in OpenAIRE",
            },
        },
        "Accessible": {
            "A1_data_retrievable_by_identifier": {
                "status": "complete",
                "evidence": "GitHub repo is public (code); large data files excluded but documented",
                "action": "None",
            },
            "A1_1_protocol_open_free_universal": {
                "status": "complete",
                "evidence": "HTTPS protocol; no authentication required for code",
                "action": "None",
            },
            "A1_2_protocol_allows_authentication_authorization": {
                "status": "complete",
                "evidence": "GitHub supports OAuth; repository access controlled",
                "action": "None",
            },
            "A2_metadata_accessible_even_data_no_longer_available": {
                "status": "partial",
                "evidence": "Metadata in JSON files committed to repo; Zenodo deposit will provide persistence",
                "action": "Complete Zenodo deposit",
            },
        },
        "Interoperable": {
            "I1_data_uses_formal_accessible_shared_language": {
                "status": "complete",
                "evidence": "CSV, JSON, and Markdown formats; no proprietary formats",
                "action": "None",
            },
            "I1_1_data_uses_fair_vocabularies": {
                "status": "partial",
                "evidence": "Gene symbols use HGNC; no formal ontology mappings (e.g., NCBITaxon, BFO)",
                "action": "Add ontology annotations: NCBITaxon for species, OBI for assays",
            },
            "I2_data_includes_qualified_references": {
                "status": "partial",
                "evidence": "References to GEO accessions in metadata; not all references qualified with persistent IDs",
                "action": "Add CrossRef DOIs for all cited publications; ORCID for authors",
            },
        },
        "Reusable": {
            "R1_metadata_richly_described": {
                "status": "complete",
                "evidence": "Comprehensive README, model card, data dictionary in output/",
                "action": "None",
            },
            "R1_1_data_released_with_clear_accessible_data_usage_license": {
                "status": "complete",
                "evidence": "MIT License in LICENSE file; .zenodo.json specifies license",
                "action": "None",
            },
            "R1_2_data_associated_with_detailed_provenance": {
                "status": "partial",
                "evidence": "Scripts document data sources; full provenance chain (W3C PROV) not yet implemented",
                "action": "Add W3C PROV metadata or use Dataverse/DataONE provenance tools",
            },
            "R1_3_data_meets_domain_relevant_community_standards": {
                "status": "partial",
                "evidence": "Follows GEO/MIAME for transcriptomics; does not yet follow MIABIS or ISA-Tab",
                "action": "Add ISA-Tab metadata for multi-omic experiments",
            },
        },
    },
    "overall_score": {
        "complete": 9,
        "partial": 7,
        "missing": 0,
        "total": 16,
        "percentage": 9 / 16 * 100 + (7 / 16 * 50),  # partial = 50% credit
    },
    "recommendations": [
        "1. Deposit code + metadata in Zenodo to obtain DOI (F1, F4, A2)",
        "2. Add HGNC and NCBITaxon ontology mappings to gene and species metadata (I1.1)",
        "3. Add ORCID identifiers for all authors in .zenodo.json (I2)",
        "4. Generate W3C PROV provenance trace for key result files (R1.2)",
        "5. Add ISA-Tab format for experimental metadata (R1.3)",
        "6. Register dataset in OpenAIRE and Google Dataset Search (F4)",
        "7. Add CITATION.cff file for standardized citation metadata (F2)",
    ],
}

# Compute percentage properly
checklist["overall_score"]["percentage"] = round(
    (checklist["overall_score"]["complete"] * 100 + checklist["overall_score"]["partial"] * 50) / checklist["overall_score"]["total"], 1
)

print("\n─── FAIR Assessment ───")
print(f"  Complete:   {checklist['overall_score']['complete']}/{checklist['overall_score']['total']}")
print(f"  Partial:    {checklist['overall_score']['partial']}/{checklist['overall_score']['total']}")
print(f"  Score:      {checklist['overall_score']['percentage']:.1f}%")

for principle, items in checklist["principles"].items():
    n_complete = sum(1 for v in items.values() if v["status"] == "complete")
    n_total = len(items)
    print(f"\n  {principle}: {n_complete}/{n_total}")
    for item, detail in items.items():
        icon = "✅" if detail["status"] == "complete" else "🟡" if detail["status"] == "partial" else "❌"
        print(f"    {icon} {item}")

print("\n─── Recommendations ───")
for r in checklist["recommendations"]:
    print(f"  {r}")

json.dump(checklist, open(OUT / "fair_compliance.json", "w"), indent=2)
print("\nSaved to fair_compliance.json")
print("DONE")
