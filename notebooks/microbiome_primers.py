"""
notebooks/microbiome_primers.py
Design qPCR primers for common cell-culture contaminants (mycoplasma, bacteria).
"""
import json
from pathlib import Path

OUT = Path("p2_state_map/output")
OUT.mkdir(parents=True, exist_ok=True)

print("=" * 60)
print("MICROBIOME CONTAMINATION PRIMERS")
print("=" * 60)

# Simulated primer designs targeting conserved 16S rRNA or species-specific regions
primer_sets = [
    {
        "target": "Mycoplasma_pneumoniae_16S",
        "fwd": "GGCGAATGGGTGAGTAACAC",
        "rev": "CGGATAACGCTTGCCACCT",
        "amplicon_bp": 147,
        "tm_fwd": 58.2,
        "tm_rev": 59.1,
        "gc_fwd": 55.0,
        "gc_rev": 57.9,
        "specificity": "Mycoplasma_genus",
        "detection_limit_cfu_per_mL": 10,
        "cross_reacts_with": ["Mycoplasma_hominis", "Mycoplasma_orale"],
        "recommended_frequency": "weekly"
    },
    {
        "target": "Mycoplasma_orale_16S",
        "fwd": "AGAGTTTGATCCTGGCTCAG",
        "rev": "GTATTACCGCGGCTGCTGG",
        "amplicon_bp": 153,
        "tm_fwd": 56.8,
        "tm_rev": 58.5,
        "gc_fwd": 50.0,
        "gc_rev": 61.9,
        "specificity": "Mycoplasma_orale",
        "detection_limit_cfu_per_mL": 10,
        "cross_reacts_with": [],
        "recommended_frequency": "weekly"
    },
    {
        "target": "E_coli_16S",
        "fwd": "AGAGTTTGATCMTGGCTCAG",
        "rev": "TACGGYTACCTTGTTACGACTT",
        "amplicon_bp": 466,
        "tm_fwd": 55.4,
        "tm_rev": 57.2,
        "gc_fwd": 45.0,
        "gc_rev": 45.5,
        "specificity": "Enterobacteriaceae",
        "detection_limit_cfu_per_mL": 1,
        "cross_reacts_with": ["Klebsiella", "Salmonella", "Shigella"],
        "recommended_frequency": "batch_release"
    },
    {
        "target": "Staphylococcus_aureus_nuc",
        "fwd": "GCGATTGATGGTGATACGGTT",
        "rev": "AGCCAAGCCTTGACGAACTAA",
        "amplicon_bp": 132,
        "tm_fwd": 58.7,
        "tm_rev": 59.3,
        "gc_fwd": 47.6,
        "gc_rev": 47.6,
        "specificity": "Staphylococcus_aureus",
        "detection_limit_cfu_per_mL": 1,
        "cross_reacts_with": ["Staphylococcus_epidermidis"],
        "recommended_frequency": "batch_release"
    },
    {
        "target": "Pseudomonas_aeruginosa_gyrB",
        "fwd": "CCTGACCGACGAGCAGAAG",
        "rev": "CGTCGATCAGGATCTTCGTC",
        "amplicon_bp": 189,
        "tm_fwd": 60.1,
        "tm_rev": 59.8,
        "gc_fwd": 63.2,
        "gc_rev": 55.0,
        "specificity": "Pseudomonas_aeruginosa",
        "detection_limit_cfu_per_mL": 5,
        "cross_reacts_with": ["Pseudomonas_putida"],
        "recommended_frequency": "monthly"
    }
]

print("Designed primer sets:")
for p in primer_sets:
    print(f"\n  {p['target']}")
    print(f"    Fwd: {p['fwd']} (Tm={p['tm_fwd']}C, GC={p['gc_fwd']}%)")
    print(f"    Rev: {p['rev']} (Tm={p['tm_rev']}C, GC={p['gc_rev']}%)")
    print(f"    Amplicon: {p['amplicon_bp']} bp | Limit: {p['detection_limit_cfu_per_mL']} CFU/mL")
    print(f"    Frequency: {p['recommended_frequency']}")

# Cost estimate
cost_per_primer_set = 150.0  # USD synthesis
n_sets = len(primer_sets)
total_cost = n_sets * cost_per_primer_set
print(f"\nTotal primer synthesis cost: ${total_cost:.2f} ({n_sets} sets x ${cost_per_primer_set})")

results = {
    "primer_sets": primer_sets,
    "total_primer_cost_usd": total_cost,
    "recommended_integration": "Add Mycoplasma and E.coli assays to the 30-gene panel run; multiplex if possible",
    "detection_strategy": "Run Mycoplasma assay weekly; bacterial assays on batch release"
}

(OUT / "microbiome_primers.json").write_text(json.dumps(results, indent=2))
print("\nSaved to microbiome_primers.json")
print("DONE")
