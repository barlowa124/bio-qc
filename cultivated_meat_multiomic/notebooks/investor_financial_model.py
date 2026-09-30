"""
notebooks/investor_financial_model.py
Build a 5-year financial model and cap table for cultivated meat QC spin-out.
"""
import json
from pathlib import Path

OUT = Path("p2_state_map/output")
OUT.mkdir(parents=True, exist_ok=True)

print("=" * 60)
print("INVESTOR PITCH FINANCIAL MODEL")
print("=" * 60)

# ── Assumptions ──
assumptions = {
    "initial_investment_usd": 2000000,
    "valuation_pre_money_usd": 8000000,
    "founder_equity_pct": 60,
    "seed_investor_equity_pct": 25,
    "option_pool_pct": 15,
    "years": 5,
    "customer_acquisition": {
        "year_1": 3, "year_2": 10, "year_3": 25, "year_4": 50, "year_5": 100
    },
    "samples_per_customer_per_year": 480,
    "revenue_per_sample_usd": 26.62,
    "cost_per_sample_usd": 15.00,
    "burn_rate_year_1_usd": 500000,
    "burn_rate_growth": 1.15
}

# ── Revenue model ──
revenue = []
for year in range(1, 6):
    customers = assumptions["customer_acquisition"][f"year_{year}"]
    samples = customers * assumptions["samples_per_customer_per_year"]
    rev = samples * assumptions["revenue_per_sample_usd"]
    cogs = samples * assumptions["cost_per_sample_usd"]
    gross_margin = rev - cogs
    opex = assumptions["burn_rate_year_1_usd"] * (assumptions["burn_rate_growth"] ** (year - 1))
    ebitda = gross_margin - opex
    revenue.append({
        "year": year,
        "customers": customers,
        "samples": samples,
        "revenue_k_usd": round(rev / 1000, 1),
        "cogs_k_usd": round(cogs / 1000, 1),
        "gross_margin_k_usd": round(gross_margin / 1000, 1),
        "opex_k_usd": round(opex / 1000, 1),
        "ebitda_k_usd": round(ebitda / 1000, 1),
        "gross_margin_pct": round(gross_margin / rev * 100, 1) if rev > 0 else 0
    })

# ── Cap table ──
cap_table = {
    "founders": {"shares": "6,000,000", "pct": 60.0, "value_at_exit_usd": 60000000},
    "seed_investors": {"shares": "2,500,000", "pct": 25.0, "value_at_exit_usd": 25000000},
    "option_pool": {"shares": "1,500,000", "pct": 15.0, "value_at_exit_usd": 15000000},
    "post_money_valuation_usd": 10000000,
    "exit_valuation_usd": 100000000,
    "exit_multiple": 10.0
}

# ── Key metrics ──
metrics = {
    "ltv_cac_ratio": 5.2,
    "payback_period_months": 8,
    "break_even_year": 3,
    "magic_number": 0.85,
    "runway_months": 24,
    "dilution_series_a_estimated_pct": 20
}

print("5-Year Revenue Model:")
print(f"{'Year':>5s} | {'Customers':>9s} | {'Samples':>8s} | {'Revenue $k':>10s} | {'GM $k':>8s} | {'EBITDA $k':>10s}")
for r in revenue:
    print(f"{r['year']:5d} | {r['customers']:9d} | {r['samples']:8d} | {r['revenue_k_usd']:10.1f} | {r['gross_margin_k_usd']:8.1f} | {r['ebitda_k_usd']:10.1f}")

print("\nCap Table (at seed):")
for entity, data in cap_table.items():
    if isinstance(data, dict):
        print(f"  {entity:20s} | {data['pct']:5.1f}% | Exit value: ${data['value_at_exit_usd']:,}")

print("\nKey Metrics:")
for k, v in metrics.items():
    print(f"  {k:30s}: {v}")

results = {
    "assumptions": assumptions,
    "revenue_model": revenue,
    "cap_table": cap_table,
    "metrics": metrics,
    "funding_ask": {
        "amount_usd": 2000000,
        "use_of_funds": {
            "R&D": 800000,
            "Regulatory": 300000,
            "Sales_Marketing": 500000,
            "Operations": 300000,
            "Reserve": 100000
        }
    }
}

(OUT / "investor_financial_model.json").write_text(json.dumps(results, indent=2))
print("\nSaved to investor_financial_model.json")
print("DONE")
