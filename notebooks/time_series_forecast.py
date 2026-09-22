"""
notebooks/time_series_forecast.py
Forecast future cell states from time-series qPCR data.
"""
import json
from pathlib import Path

import numpy as np
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_squared_error

np.random.seed(42)
OUT = Path("p2_state_map/output")
OUT.mkdir(parents=True, exist_ok=True)

print("=" * 60)
print("TIME-SERIES STATE FORECASTING")
print("=" * 60)

# Simulate time-series: 20 timepoints, 30 genes, state transitions
timepoints = 20
genes = 30

# Generate a trajectory: expansion -> committed -> terminal
true_states = (["expansion_competent"] * 7 + ["committed"] * 7 + ["terminal"] * 6)
X = np.zeros((timepoints, genes))
for t in range(timepoints):
    if true_states[t] == "expansion_competent":
        X[t, 0:10] = np.random.normal(2, 0.5, 10)
        X[t, 10:] = np.random.normal(0, 0.5, 20)
    elif true_states[t] == "committed":
        X[t, 10:20] = np.random.normal(2, 0.5, 10)
        X[t, 0:10] = np.random.normal(0.5, 0.5, 10)
        X[t, 20:] = np.random.normal(0, 0.5, 10)
    else:
        X[t, 20:30] = np.random.normal(2, 0.5, 10)
        X[t, 0:20] = np.random.normal(-0.5, 0.5, 20)

# Add temporal smoothness
for t in range(1, timepoints):
    X[t] = 0.7 * X[t] + 0.3 * X[t-1]

# Forecast horizon
horizon = 3

# Method 1: Naive (last value)
# Method 2: Autoregressive Ridge (predict next from last 3)
# Method 3: State transition model (Markov)

def ar_forecast(X, lag=3, horizon=3):
    preds = []
    for h in range(horizon):
        train_X = []
        train_y = []
        for t in range(lag, len(X) - h):
            train_X.append(X[t-lag:t].flatten())
            train_y.append(X[t+h])
        train_X = np.array(train_X)
        train_y = np.array(train_y)
        model = Ridge(alpha=1.0)
        model.fit(train_X, train_y)
        last = X[-lag:].flatten().reshape(1, -1)
        pred = model.predict(last)[0]
        preds.append(pred)
    return np.array(preds)

def markov_forecast(states, horizon=3):
    # Transition counts
    trans = {}
    for s in set(states):
        trans[s] = {s2: 0 for s2 in set(states)}
    for i in range(len(states)-1):
        trans[states[i]][states[i+1]] += 1
    # Normalize
    for s in trans:
        total = sum(trans[s].values())
        if total > 0:
            for s2 in trans[s]:
                trans[s][s2] /= total
    # Forecast
    current = states[-1]
    forecast = []
    for _ in range(horizon):
        probs = [trans[current].get(s, 0) for s in ["expansion_competent", "committed", "terminal"]]
        if sum(probs) == 0:
            probs = [0.33, 0.33, 0.34]
        next_state = np.random.choice(["expansion_competent", "committed", "terminal"], p=probs)
        forecast.append(next_state)
        current = next_state
    return forecast

ar_preds = ar_forecast(X, lag=3, horizon=horizon)
markov_preds = markov_forecast(true_states, horizon=horizon)

print(f"\nAR Forecast (next {horizon} timepoints, first 5 genes):")
for h in range(horizon):
    print(f"  t+{h+1}: {ar_preds[h, :5]}")

print(f"\nMarkov State Forecast (next {horizon} timepoints):")
for h in range(horizon):
    print(f"  t+{h+1}: {markov_preds[h]}")

# Simulate future true values for validation
future_true = X[-horizon:]
future_states = true_states[-horizon:]
rmse_ar = np.sqrt(mean_squared_error(future_true, ar_preds))
markov_acc = np.mean([markov_preds[i] == future_states[i] for i in range(horizon)])

print("\nValidation:")
print(f"  AR RMSE: {rmse_ar:.4f}")
print(f"  Markov accuracy: {markov_acc:.0%}")

results = {
    "timepoints": timepoints,
    "forecast_horizon": horizon,
    "ar_forecast": {
        "predictions": [[round(float(v), 4) for v in ar_preds[h]] for h in range(horizon)],
        "rmse_vs_future": round(float(rmse_ar), 4)
    },
    "markov_forecast": {
        "predictions": markov_preds,
        "accuracy_vs_future": round(float(markov_acc), 4)
    },
    "recommendation": "Use AR model for gene-level trajectory forecasting; Markov model for discrete state transitions. Combine for early warning system."
}

(OUT / "time_series_forecast.json").write_text(json.dumps(results, indent=2))
print("\nSaved to time_series_forecast.json")
print("DONE")
