"""Techno-economic analysis (TEA) for 30-gene qPCR panel at manufacturing scale."""
import json
from pathlib import Path

import numpy as np

PROJ = Path(__file__).resolve().parents[1]
OUT = PROJ / "p2_state_map/output"
OUT.mkdir(parents=True, exist_ok=True)

print("=" * 60)
print("TECHNO-ECONOMIC ANALYSIS (TEA)")
print("=" * 60)

# ── Base parameters ──
PANEL_COST_PER_BATCH = 75.0  # $50-100 midpoint
GENES = 30
SAMPLES_PER_BATCH = 8  # 96-well plate, 30 genes + controls + replicates
TURNAROUND_HOURS = 5.0  # 4-6 hour midpoint
LABOR_RATE = 45.0  # $/hour technician

# ── Scenario: lab-scale QC (current) ──
print("\n─── Scenario 1: Lab-Scale QC ───")
lab = {
    "batches_per_week": 5,
    "samples_per_batch": SAMPLES_PER_BATCH,
    "panel_cost_per_batch": PANEL_COST_PER_BATCH,
    "labor_hours_per_batch": 2.0,
    "labor_rate": LABOR_RATE,
    "equipment_depreciation_per_year": 5000.0,
    "overhead_multiplier": 1.3,
}
lab["samples_per_week"] = lab["batches_per_week"] * lab["samples_per_batch"]
lab["panel_cost_per_year"] = lab["batches_per_week"] * 52 * lab["panel_cost_per_batch"]
lab["labor_cost_per_year"] = lab["batches_per_week"] * 52 * lab["labor_hours_per_batch"] * lab["labor_rate"]
lab["equipment_per_year"] = lab["equipment_depreciation_per_year"]
lab["total_cost_per_year"] = (lab["panel_cost_per_year"] + lab["labor_cost_per_year"] + lab["equipment_per_year"]) * lab["overhead_multiplier"]
lab["cost_per_sample"] = lab["total_cost_per_year"] / (lab["samples_per_week"] * 52)
lab["cost_per_batch"] = lab["total_cost_per_year"] / (lab["batches_per_week"] * 52)
for k, v in lab.items():
    if isinstance(v, float):
        print(f"  {k}: ${v:,.2f}" if "cost" in k or "_per_" in k else f"  {k}: {v:,.1f}")

# ── Scenario: pilot plant (1000L bioreactor) ──
print("\n─── Scenario 2: Pilot Plant (1000L) ───")
pilot = {
    "bioreactor_volume_L": 1000,
    "cell_density_M_per_mL": 10.0,  # 10 million cells/mL at harvest
    "samples_per_run": 20,  # QC samples across run
    "runs_per_month": 4,
    "panel_cost_per_batch": PANEL_COST_PER_BATCH * 0.85,  # modest volume discount
    "labor_hours_per_batch": 1.5,
    "labor_rate": LABOR_RATE,
    "equipment_depreciation_per_year": 15000.0,
    "overhead_multiplier": 1.2,
}
pilot["batches_per_year"] = pilot["runs_per_month"] * 12
pilot["samples_per_year"] = pilot["batches_per_year"] * pilot["samples_per_run"]
pilot["panel_cost_per_year"] = pilot["batches_per_year"] * pilot["panel_cost_per_batch"]
pilot["labor_cost_per_year"] = pilot["batches_per_year"] * pilot["labor_hours_per_batch"] * pilot["labor_rate"]
pilot["equipment_per_year"] = pilot["equipment_depreciation_per_year"]
pilot["total_cost_per_year"] = (pilot["panel_cost_per_year"] + pilot["labor_cost_per_year"] + pilot["equipment_per_year"]) * pilot["overhead_multiplier"]
pilot["cost_per_sample"] = pilot["total_cost_per_year"] / pilot["samples_per_year"]
pilot["cost_per_1000L_run"] = pilot["total_cost_per_year"] / pilot["batches_per_year"]
pilot["cost_per_million_cells"] = pilot["total_cost_per_year"] / (pilot["bioreactor_volume_L"] * pilot["cell_density_M_per_mL"] * pilot["batches_per_year"])
for k, v in pilot.items():
    if isinstance(v, float):
        print(f"  {k}: ${v:,.4f}" if "cost" in k else f"  {k}: {v:,.2f}")

# ── Scenario: manufacturing (10,000L) ──
print("\n─── Scenario 3: Manufacturing (10,000L) ───")
mfg = {
    "bioreactor_volume_L": 10000,
    "cell_density_M_per_mL": 15.0,
    "samples_per_run": 50,
    "runs_per_month": 20,
    "panel_cost_per_batch": PANEL_COST_PER_BATCH * 0.60,  # high-volume discount
    "labor_hours_per_batch": 1.0,  # automation
    "labor_rate": LABOR_RATE * 1.5,  # senior technician
    "equipment_depreciation_per_year": 50000.0,
    "overhead_multiplier": 1.15,
    "automation_capital": 200000.0,  # liquid handler + qPCR robot
    "automation_amortization_years": 5,
}
mfg["batches_per_year"] = mfg["runs_per_month"] * 12
mfg["samples_per_year"] = mfg["batches_per_year"] * mfg["samples_per_run"]
mfg["panel_cost_per_year"] = mfg["batches_per_year"] * mfg["panel_cost_per_batch"]
mfg["labor_cost_per_year"] = mfg["batches_per_year"] * mfg["labor_hours_per_batch"] * mfg["labor_rate"]
mfg["equipment_per_year"] = mfg["equipment_depreciation_per_year"] + (mfg["automation_capital"] / mfg["automation_amortization_years"])
mfg["total_cost_per_year"] = (mfg["panel_cost_per_year"] + mfg["labor_cost_per_year"] + mfg["equipment_per_year"]) * mfg["overhead_multiplier"]
mfg["cost_per_sample"] = mfg["total_cost_per_year"] / mfg["samples_per_year"]
mfg["cost_per_10000L_run"] = mfg["total_cost_per_year"] / mfg["batches_per_year"]
mfg["cost_per_million_cells"] = mfg["total_cost_per_year"] / (mfg["bioreactor_volume_L"] * mfg["cell_density_M_per_mL"] * mfg["batches_per_year"])
mfg["cost_per_kg_meat"] = mfg["total_cost_per_year"] / (mfg["bioreactor_volume_L"] * mfg["cell_density_M_per_mL"] * mfg["batches_per_year"] / 1000.0 * 0.3)  # ~30% yield to meat
for k, v in mfg.items():
    if isinstance(v, float):
        print(f"  {k}: ${v:,.4f}" if "cost" in k else f"  {k}: {v:,.2f}")

# ── Sensitivity: what if panel cost drops? ──
print("\n─── Sensitivity: Panel Cost vs Total Cost ───")
sensitivity = []
for cost in [30, 40, 50, 60, 75, 100]:
    panel_year = mfg["batches_per_year"] * cost
    total = (panel_year + mfg["labor_cost_per_year"] + mfg["equipment_per_year"]) * mfg["overhead_multiplier"]
    per_sample = total / mfg["samples_per_year"]
    sensitivity.append({"panel_cost_per_batch": cost, "total_per_year": total, "cost_per_sample": per_sample})
    print(f"  ${cost}/batch → ${per_sample:.2f}/sample")

# ── Comparison with alternative QC methods ──
print("\n─── Alternative QC Methods ───")
alternatives = {
    "30_gene_qPCR_panel": {"cost_per_sample": pilot["cost_per_sample"], "turnaround_hours": 5, "accuracy": 0.967, "throughput": "high"},
    "bulk_RNA_seq": {"cost_per_sample": 250.0, "turnaround_hours": 336, "accuracy": 0.98, "throughput": "medium"},
    "snRNA_seq": {"cost_per_sample": 800.0, "turnaround_hours": 504, "accuracy": 0.95, "throughput": "low"},
    "proteomics": {"cost_per_sample": 400.0, "turnaround_hours": 168, "accuracy": 0.92, "throughput": "medium"},
    "metabolomics": {"cost_per_sample": 300.0, "turnaround_hours": 120, "accuracy": 0.88, "throughput": "medium"},
    "microscopy_morphology": {"cost_per_sample": 15.0, "turnaround_hours": 2, "accuracy": 0.75, "throughput": "very_high"},
    "ATP_luminescence": {"cost_per_sample": 5.0, "turnaround_hours": 1, "accuracy": 0.70, "throughput": "very_high"},
}
for name, data in alternatives.items():
    print(f"  {name}: ${data['cost_per_sample']:.2f}/sample, {data['turnaround_hours']}h, {data['accuracy']:.0%} accuracy")

# ── Save ──
results = {
    "lab_scale": {k: float(v) if isinstance(v, np.floating) else v for k, v in lab.items()},
    "pilot_plant": {k: float(v) if isinstance(v, np.floating) else v for k, v in pilot.items()},
    "manufacturing": {k: float(v) if isinstance(v, np.floating) else v for k, v in mfg.items()},
    "sensitivity": sensitivity,
    "alternatives": alternatives,
    "panel_genes": 30,
    "genes": ["UXS1", "PLOD1", "MALAT1", "C1D", "KIF1B", "UPK1B", "SNHG3", "C1QBP", "LMNA", "TOMM7",
              "MRPL32", "AEBP1", "PLXNA1", "GLG1", "CTSA", "NPC2", "RAB42", "C3orf72", "ZNF527", "FCRLA",
              "C11orf63", "EMC1", "HRH4", "PLEKHG4B", "PPIEL", "XKR9", "APOL1", "LINC00574", "UGT8", "IFITM3"],
}
json.dump(results, open(OUT / "techno_economic_analysis.json", "w"), indent=2)
print("\nSaved to techno_economic_analysis.json")
print("DONE")
