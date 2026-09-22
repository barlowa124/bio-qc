"""
notebooks/supply_chain_risk.py
Map primer/reagent availability, geographic sourcing, and lead-time risks.
"""
import json
from pathlib import Path

OUT = Path("p2_state_map/output")
OUT.mkdir(parents=True, exist_ok=True)

print("=" * 60)
print("SUPPLY CHAIN RISK ANALYSIS")
print("=" * 60)

# ── Primer supplier landscape ──
suppliers = [
    {"name": "IDT (Integrated DNA Technologies)", "country": "USA", "lead_time_days": 5, "reliability": 0.99, "share": 0.35},
    {"name": "Thermo Fisher (Invitrogen)", "country": "USA", "lead_time_days": 7, "reliability": 0.98, "share": 0.30},
    {"name": "Sigma-Aldrich (Merck)", "country": "Germany", "lead_time_days": 10, "reliability": 0.97, "share": 0.20},
    {"name": "Eurofins Genomics", "country": "Germany", "lead_time_days": 8, "reliability": 0.96, "share": 0.10},
    {"name": "GenScript", "country": "China", "lead_time_days": 12, "reliability": 0.94, "share": 0.05},
]

# ── Reagent categories ──
reagents = [
    {"category": "TaqMan probes / primers", "critical": True, "suppliers": ["IDT", "Thermo", "Sigma"], "alt_available": True},
    {"category": "RT-qPCR master mix", "critical": True, "suppliers": ["Thermo", "Roche", "Bio-Rad"], "alt_available": True},
    {"category": "RNA extraction kits", "critical": True, "suppliers": ["Qiagen", "Thermo", "Zymo"], "alt_available": True},
    {"category": "Cell culture media (DMEM/F12)", "critical": True, "suppliers": ["Thermo", "Sigma", "Corning"], "alt_available": True},
    {"category": "Fetal bovine serum", "critical": True, "suppliers": ["Thermo", "Sigma", "HyClone"], "alt_available": False},
    {"category": "Bioreactor consumables", "critical": True, "suppliers": ["Sartorius", "Thermo", "Eppendorf"], "alt_available": False},
    {"category": "Sequencing reagents (if backup RNA-seq)", "critical": False, "suppliers": ["Illumina", "Thermo"], "alt_available": False},
]

# ── Risk scenarios ──
scenarios = [
    {"name": "Base case", "probability": 1.0, "impact_days": 0},
    {"name": "Single supplier delay (IDT)", "probability": 0.15, "impact_days": 5},
    {"name": "Single supplier delay (Thermo)", "probability": 0.15, "impact_days": 7},
    {"name": "Pandemic / border closure", "probability": 0.05, "impact_days": 30},
    {"name": "Trade tariff (China-USA)", "probability": 0.10, "impact_days": 14, "cost_increase_pct": 15},
    {"name": "Regional conflict (EU)", "probability": 0.03, "impact_days": 21},
]

# ── Mitigation strategies ──
mitigations = [
    "Dual-source all critical primers (IDT + Thermo)",
    "Maintain 90-day safety stock of master mix and extraction kits",
    "Pre-validate alternative suppliers with identical primer sequences",
    "Regional diversification: keep 30% inventory in Singapore for APAC manufacturing",
    "Long-term contracts with volume commitments for price stability",
    "Investigate in-house primer synthesis for highest-turnover assays"
]

# ── Calculations ──
total_panel_cost = 266.20  # 30 genes x ~$8.87 avg (from TEA)
avg_lead = sum(s["lead_time_days"] * s["share"] for s in suppliers)
concentration_risk = sum(s["share"]**2 for s in suppliers)

print("Supplier landscape:")
for s in suppliers:
    print(f"  {s['name']:35s} | {s['country']:10s} | Lead: {s['lead_time_days']}d | Reliability: {s['reliability']:.0%} | Share: {s['share']:.0%}")

print(f"\nWeighted average lead time: {avg_lead:.1f} days")
print(f"Herfindahl concentration index: {concentration_risk:.2f} (0=diverse, 1=monopoly)")

print("\nReagent risk matrix:")
for r in reagents:
    risk = "HIGH" if r["critical"] and not r["alt_available"] else "MEDIUM" if r["critical"] else "LOW"
    print(f"  {r['category']:35s} | Critical={r['critical']} | Alt={r['alt_available']} | Risk={risk}")

print("\nMitigation strategies:")
for i, m in enumerate(mitigations, 1):
    print(f"  {i}. {m}")

results = {
    "suppliers": suppliers,
    "reagents": reagents,
    "scenarios": scenarios,
    "mitigations": mitigations,
    "metrics": {
        "weighted_avg_lead_time_days": round(avg_lead, 1),
        "concentration_index": round(concentration_risk, 3),
        "total_panel_cost_usd": total_panel_cost
    },
    "recommendation": "Dual-source all critical primers; maintain 90-day safety stock; regional inventory in Singapore"
}

(OUT / "supply_chain_risk.json").write_text(json.dumps(results, indent=2))
print("\nSaved to supply_chain_risk.json")
print("DONE")
