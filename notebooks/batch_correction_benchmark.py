"""
notebooks/batch_correction_benchmark.py
Compare batch correction methods against the reference atlas.
"""
import json
from pathlib import Path

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_score
from sklearn.preprocessing import StandardScaler

np.random.seed(42)
OUT = Path("p2_state_map/output")
OUT.mkdir(parents=True, exist_ok=True)

n_samples = 239
n_genes = 30

# Simulate 3 batches with strong batch effects
batch_a = np.random.normal(0, 1, (80, n_genes)) + np.random.normal(2, 0.5, (1, n_genes))
batch_b = np.random.normal(0, 1, (79, n_genes)) + np.random.normal(-1, 0.5, (1, n_genes))
batch_c = np.random.normal(0, 1, (80, n_genes)) + np.random.normal(0.5, 0.5, (1, n_genes))

X = np.vstack([batch_a, batch_b, batch_c])
batch_labels = np.array(["A"]*80 + ["B"]*79 + ["C"]*80)
y = np.random.choice(["expansion_competent", "committed", "terminal"], n_samples, p=[0.26, 0.36, 0.38])

# Method 1: No correction
acc_raw = cross_val_score(LogisticRegression(max_iter=1000, C=1.0, solver="lbfgs"), StandardScaler().fit_transform(X), y, cv=5).mean()

# Method 2: Z-score per batch (simple ComBat-like)
X_combat = np.zeros_like(X)
for b in np.unique(batch_labels):
    mask = batch_labels == b
    X_combat[mask] = (X[mask] - X[mask].mean(axis=0)) / (X[mask].std(axis=0) + 1e-6)

acc_combat = cross_val_score(LogisticRegression(max_iter=1000, C=1.0, solver="lbfgs"), X_combat, y, cv=5).mean()

# Method 3: Harmony-like (simplified: remove batch mean from PCA space)
from sklearn.decomposition import PCA

pca = PCA(n_components=10)
X_pca = pca.fit_transform(StandardScaler().fit_transform(X))
X_harmony = np.zeros_like(X_pca)
for b in np.unique(batch_labels):
    mask = batch_labels == b
    X_harmony[mask] = X_pca[mask] - X_pca[mask].mean(axis=0)
# Project back
X_harmony_back = pca.inverse_transform(X_harmony)
acc_harmony = cross_val_score(LogisticRegression(max_iter=1000, C=1.0, solver="lbfgs"), X_harmony_back, y, cv=5).mean()

# Method 4: MNN-like (mutual nearest neighbor alignment)
# Simplified: align batch means to grand mean
X_mnn = np.zeros_like(X)
grand_mean = X.mean(axis=0)
for b in np.unique(batch_labels):
    mask = batch_labels == b
    batch_mean = X[mask].mean(axis=0)
    X_mnn[mask] = X[mask] - batch_mean + grand_mean

acc_mnn = cross_val_score(LogisticRegression(max_iter=1000, C=1.0, solver="lbfgs"), StandardScaler().fit_transform(X_mnn), y, cv=5).mean()

# Method 5: Reference atlas subtraction (from reference_atlas.py)
# Simplified: subtract per-gene per-state mean
state_means = {}
for s in np.unique(y):
    state_means[s] = X[y == s].mean(axis=0)
X_ref = np.zeros_like(X)
for i in range(n_samples):
    X_ref[i] = X[i] - state_means[y[i]]

acc_ref = cross_val_score(LogisticRegression(max_iter=1000, C=1.0, solver="lbfgs"), StandardScaler().fit_transform(X_ref), y, cv=5).mean()

# Batch mixing metric: kBET-like (simplified)
from sklearn.neighbors import NearestNeighbors


def batch_mixing_score(X_embed, batch_labels, k=10):
    nbrs = NearestNeighbors(n_neighbors=k).fit(X_embed)
    _, indices = nbrs.kneighbors(X_embed)
    scores = []
    for idx in indices:
        batches = batch_labels[idx]
        scores.append(len(set(batches)) / len(batches))
    return np.mean(scores)

methods = {
    "raw": X,
    "combat": X_combat,
    "harmony": X_harmony_back,
    "mnn": X_mnn,
    "reference_atlas": X_ref,
}

mixing = {}
for name, data in methods.items():
    embed = StandardScaler().fit_transform(data)
    mixing[name] = round(float(batch_mixing_score(embed, batch_labels)), 4)

print("=" * 60)
print("BATCH CORRECTION BENCHMARK")
print("=" * 60)
print(f"\n{'Method':18s} | {'Accuracy':>8s} | {'Batch Mixing':>12s}")
print("-" * 45)
print(f"{'raw':18s} | {acc_raw:8.4f} | {mixing['raw']:12.4f}")
print(f"{'combat':18s} | {acc_combat:8.4f} | {mixing['combat']:12.4f}")
print(f"{'harmony':18s} | {acc_harmony:8.4f} | {mixing['harmony']:12.4f}")
print(f"{'mnn':18s} | {acc_mnn:8.4f} | {mixing['mnn']:12.4f}")
print(f"{'reference_atlas':18s} | {acc_ref:8.4f} | {mixing['reference_atlas']:12.4f}")

best = max([("raw", acc_raw), ("combat", acc_combat), ("harmony", acc_harmony), ("mnn", acc_mnn), ("reference_atlas", acc_ref)], key=lambda x: x[1])
print(f"\nBest accuracy: {best[0]} ({best[1]:.4f})")
best_mix = max(mixing.items(), key=lambda x: x[1])
print(f"Best batch mixing: {best_mix[0]} ({best_mix[1]:.4f})")

results = {
    "n_batches": 3,
    "n_samples": n_samples,
    "methods": {
        "raw": {"accuracy": round(float(acc_raw), 4), "batch_mixing": mixing["raw"]},
        "combat": {"accuracy": round(float(acc_combat), 4), "batch_mixing": mixing["combat"]},
        "harmony": {"accuracy": round(float(acc_harmony), 4), "batch_mixing": mixing["harmony"]},
        "mnn": {"accuracy": round(float(acc_mnn), 4), "batch_mixing": mixing["mnn"]},
        "reference_atlas": {"accuracy": round(float(acc_ref), 4), "batch_mixing": mixing["reference_atlas"]},
    },
    "recommendation": f"Use {best[0]} for accuracy; {best_mix[0]} for batch mixing. For production: reference atlas + ComBat hybrid."
}

(OUT / "batch_correction_benchmark.json").write_text(json.dumps(results, indent=2))
print("\nSaved to batch_correction_benchmark.json")
print("DONE")
