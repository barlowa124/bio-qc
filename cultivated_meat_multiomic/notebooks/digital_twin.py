"""
notebooks/digital_twin.py
Bioreactor digital twin with 30-gene QC panel as closed-loop feedback controller.
"""
import json
import random
from pathlib import Path

import numpy as np

random.seed(42)
np.random.seed(42)

OUT = Path("p2_state_map/output")
OUT.mkdir(parents=True, exist_ok=True)

# ── Bioreactor ODE parameters ──
mu_max = 0.035          # 1/hr
K_s = 2.5               # g/L glucose
Y_xs = 0.45             # g cells / g glucose
Y_ws = 0.15             # g lactate / g glucose
k_d = 0.003             # 1/hr (death rate)
V = 1000                # L (pilot scale)
feed_glucose = 40       # g/L feed concentration
D_base = 0.01           # 1/hr base dilution rate

# ── QC panel parameters ──
sample_interval = 24    # hours between qPCR samples
panel_turnaround = 5    # hours from sample to result
genes = 30
cost_per_sample = 26.62 # USD from TEA

# ── State map: 0=expansion_competent, 1=committed, 2=terminal ──
STATES = ["expansion_competent", "committed", "terminal"]

def logistic_state(X, N):
    """Simple heuristic: low density + high nutrient = expansion; high density = terminal"""
    if X < 0.5e6 and N > 10:
        return 0
    elif X > 3.0e6:
        return 2
    else:
        return 1

def bioreactor_step(X, N, W, D, dt=1.0):
    """Euler step for chemostat-like bioreactor"""
    mu = mu_max * N / (K_s + N)
    dXdt = (mu - k_d) * X - D * X
    dNdt = - (mu * X) / Y_xs + D * (feed_glucose - N)
    dWdt = (mu * X) / Y_xs * Y_ws - D * W
    X = max(0.0, X + dXdt * dt)
    N = max(0.0, N + dNdt * dt)
    W = max(0.0, W + dWdt * dt)
    return X, N, W

# ── Closed-loop simulation ──
print("=" * 60)
print("DIGITAL TWIN: BIOREACTOR + QC PANEL FEEDBACK")
print("=" * 60)

results = []
for scenario_name, D0 in [("batch", 0.0), ("fed_batch", 0.005), ("perfusion", 0.02)]:
    X, N, W = 0.1e6, 40.0, 0.1  # initial cells/mL, glucose g/L, lactate g/L
    t = 0
    history = []
    total_cost = 0
    interventions = 0
    harvested = False
    D = D0

    while t < 336 and not harvested:  # 14 days max
        X, N, W = bioreactor_step(X, N, W, D)
        t += 1

        # QC panel sampling
        if t % sample_interval == 0:
            true_state = logistic_state(X, N)
            # Add measurement noise
            measured_state = true_state
            if np.random.rand() < 0.05:
                measured_state = random.choice([s for s in [0,1,2] if s != true_state])
            total_cost += cost_per_sample

            # Feedback controller
            action = "maintain"
            if measured_state == 2:  # terminal
                action = "harvest"
                harvested = True
            elif measured_state == 0 and D < 0.015:  # expansion competent, increase nutrient
                D = min(0.03, D + 0.002)
                action = "increase_feed"
                interventions += 1
            elif measured_state == 1 and W > 8.0:  # committed + high waste, dilute
                D = min(0.03, D + 0.003)
                action = "dilute_waste"
                interventions += 1

            history.append({
                "hour": int(t),
                "cell_density_M_per_mL": round(X / 1e6, 3),
                "glucose_g_per_L": round(N, 3),
                "lactate_g_per_L": round(W, 3),
                "dilution_rate_1_per_hr": round(D, 5),
                "true_state": STATES[true_state],
                "measured_state": STATES[measured_state],
                "controller_action": action,
                "cumulative_cost": round(total_cost, 2)
            })

    results.append({
        "scenario": scenario_name,
        "final_density_M_per_mL": round(X / 1e6, 3),
        "final_glucose_g_per_L": round(N, 3),
        "final_lactate_g_per_L": round(W, 3),
        "harvested": harvested,
        "harvest_time_hr": t if harvested else None,
        "total_qc_cost_usd": round(total_cost, 2),
        "interventions": interventions,
        "history": history
    })
    print(f"\nScenario: {scenario_name}")
    print(f"  Harvested: {harvested} at t={t}hr" if harvested else f"  Not harvested by t={t}hr")
    print(f"  Final density: {X/1e6:.2f} M cells/mL")
    print(f"  Total QC cost: ${total_cost:.2f}")
    print(f"  Interventions: {interventions}")

# ── Summary metrics ──
print("\n" + "=" * 60)
print("SUMMARY")
print("=" * 60)
for r in results:
    print(f"{r['scenario']:12s} | Harvest: {str(r['harvested']):5s} | Cost: ${r['total_qc_cost_usd']:.2f} | Interventions: {r['interventions']}")

(OUT / "digital_twin.json").write_text(json.dumps(results, indent=2))
print("\nSaved to digital_twin.json")
print("DONE")
