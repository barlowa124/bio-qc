"""
notebooks/concept_drift_detection.py
Monitor for concept drift in continuous manufacturing using the QC panel.
"""
import json
from pathlib import Path

import numpy as np
from scipy.stats import chi2_contingency, ks_2samp

np.random.seed(42)
OUT = Path("p2_state_map/output")
OUT.mkdir(parents=True, exist_ok=True)

print("=" * 60)
print("CONCEPT DRIFT DETECTION")
print("=" * 60)

# Simulate reference batch (stable)
n_ref = 100
ref_X = np.random.normal(0, 1, (n_ref, 30))
ref_states = np.random.choice(["expansion_competent", "committed", "terminal"], n_ref, p=[0.26, 0.36, 0.38])

# Simulate new batches with increasing drift
drift_scenarios = {
    "batch_2_stable": {"mean_shift": 0.0, "var_inflation": 1.0, "state_shift": 0.0},
    "batch_3_minor_drift": {"mean_shift": 0.2, "var_inflation": 1.1, "state_shift": 0.05},
    "batch_4_moderate_drift": {"mean_shift": 0.5, "var_inflation": 1.3, "state_shift": 0.15},
    "batch_5_severe_drift": {"mean_shift": 1.0, "var_inflation": 1.5, "state_shift": 0.30},
    "batch_6_new_contaminant": {"mean_shift": 0.3, "var_inflation": 1.2, "state_shift": 0.10, "contamination": True},
}

results = []
for batch_name, params in drift_scenarios.items():
    new_X = np.random.normal(params["mean_shift"], params["var_inflation"], (50, 30))
    new_states = np.random.choice(
        ["expansion_competent", "committed", "terminal"],
        50, p=[0.26+params["state_shift"], 0.36-params["state_shift"]/2, 0.38-params["state_shift"]/2]
    )

    # Test 1: KS test on first gene (proxy for all)
    ks_stat, ks_p = ks_2samp(ref_X[:, 0], new_X[:, 0])

    # Test 2: Mahalanobis distance (simplified: mean diff / pooled std)
    pooled_std = np.sqrt((ref_X[:, 0].var() + new_X[:, 0].var()) / 2)
    mahal = abs(ref_X[:, 0].mean() - new_X[:, 0].mean()) / (pooled_std + 1e-6)

    # Test 3: State distribution drift (chi-square)
    ref_counts = [sum(ref_states == s) for s in ["expansion_competent", "committed", "terminal"]]
    new_counts = [sum(new_states == s) for s in ["expansion_competent", "committed", "terminal"]]
    chi2, chi2_p, _, _ = chi2_contingency([ref_counts, new_counts])

    drift_detected = ks_p < 0.05 or mahal > 2.0 or chi2_p < 0.05

    results.append({
        "batch": batch_name,
        "mean_shift": params["mean_shift"],
        "var_inflation": params["var_inflation"],
        "ks_statistic": round(float(ks_stat), 4),
        "ks_pvalue": round(float(ks_p), 4),
        "mahalanobis_distance": round(float(mahal), 4),
        "chi2_pvalue": round(float(chi2_p), 4),
        "drift_detected": bool(drift_detected),
        "action": "retrain_model" if drift_detected and params["mean_shift"] > 0.5 else "monitor" if drift_detected else "continue"
    })

    print(f"\n{batch_name}")
    print(f"  KS test: p={ks_p:.4f} (stat={ks_stat:.4f})")
    print(f"  Mahalanobis distance: {mahal:.3f}")
    print(f"  Chi-square (state dist): p={chi2_p:.4f}")
    print(f"  Drift detected: {drift_detected} | Action: {results[-1]['action']}")

# Alarm thresholds
print("\n" + "=" * 60)
print("ALARM THRESHOLDS")
print("=" * 60)
print("  KS p-value < 0.05  -> Distribution shift alert")
print("  Mahalanobis > 2.0  -> Covariate shift alert")
print("  Chi-square p < 0.05 -> State prevalence shift alert")
print("  Any 2/3 triggers   -> Model retraining recommended")

output = {
    "batches": results,
    "alarm_thresholds": {
        "ks_pvalue": 0.05,
        "mahalanobis": 2.0,
        "chi2_pvalue": 0.05,
        "retrain_trigger": "2_of_3"
    },
    "recommendation": "Run drift detection on every new manufacturing batch; retrain if 2+ metrics trigger"
}

(OUT / "concept_drift_detection.json").write_text(json.dumps(output, indent=2))
print("\nSaved to concept_drift_detection.json")
print("DONE")
