"""
notebooks/spatial_transcriptomics.py
Simulate spatial transcriptomics for muscle tissue architecture validation.
"""
import json
from pathlib import Path

import numpy as np

np.random.seed(42)
OUT = Path("p2_state_map/output")
OUT.mkdir(parents=True, exist_ok=True)

print("=" * 60)
print("SPATIAL TRANSCRIPTOMICS SIMULATION")
print("=" * 60)

# Simulate a 10x10 spatial grid representing a muscle cross-section
grid_size = 10
n_spots = grid_size * grid_size

# Each spot has proportions of 3 zones: periphery (stem), intermediate (proliferating), core (differentiated)
zone_map = np.zeros((grid_size, grid_size, 3))
for i in range(grid_size):
    for j in range(grid_size):
        dist_from_edge = min(i, j, grid_size-1-i, grid_size-1-j)
        if dist_from_edge <= 1:
            zone_map[i, j] = [0.7, 0.2, 0.1]  # periphery: stem-rich
        elif dist_from_edge <= 3:
            zone_map[i, j] = [0.2, 0.6, 0.2]  # intermediate: proliferating
        else:
            zone_map[i, j] = [0.05, 0.15, 0.8]  # core: differentiated

# Gene expression per zone (3 zones x 30 genes)
zone_expr = np.array([
    [2.0]*5 + [0.5]*10 + [0.2]*15,   # periphery: stem markers high
    [0.5]*5 + [2.0]*5 + [0.5]*5 + [0.2]*15,  # intermediate
    [0.2]*10 + [2.0]*5 + [0.5]*5 + [0.2]*10   # core: differentiation markers high
]) + np.random.normal(0, 0.3, (3, 30))

# Compute spatial expression
spatial_expr = np.zeros((n_spots, 30))
spot_coords = []
for i in range(grid_size):
    for j in range(grid_size):
        idx = i * grid_size + j
        spatial_expr[idx] = zone_map[i, j].dot(zone_expr)
        spot_coords.append((i, j))

# Identify spots where 30-gene panel would correctly classify the dominant zone
def classify_spot(expr):
    scores = [
        np.dot(expr[:5], [1]*5) / 5,   # stem score
        np.dot(expr[5:10], [1]*5) / 5,  # proliferating score
        np.dot(expr[10:15], [1]*5) / 5  # differentiated score
    ]
    return np.argmax(scores)

pred_zones = [classify_spot(spatial_expr[idx]) for idx in range(n_spots)]
true_zones = [np.argmax(zone_map[i, j]) for i in range(grid_size) for j in range(grid_size)]
accuracy = np.mean([p == t for p, t in zip(pred_zones, true_zones)])

print(f"\nSpatial classification accuracy: {accuracy:.1%}")
print(f"Total spots: {n_spots}")
print(f"Zone distribution: Periphery={sum(1 for z in true_zones if z==0)}, Intermediate={sum(1 for z in true_zones if z==1)}, Core={sum(1 for z in true_zones if z==2)}")

# Moran's I (spatial autocorrelation) for first gene
from scipy.spatial.distance import pdist, squareform

coords = np.array(spot_coords)
dist_matrix = squareform(pdist(coords))
weights = 1 / (dist_matrix + np.eye(n_spots))
weights = weights / weights.sum(axis=1, keepdims=True)

gene_vals = spatial_expr[:, 0]
gene_mean = gene_vals.mean()
moran_numer = np.sum(weights * np.outer(gene_vals - gene_mean, gene_vals - gene_mean))
moran_denom = np.sum((gene_vals - gene_mean) ** 2)
moran_i = (n_spots / weights.sum()) * (moran_numer / moran_denom)

print(f"Moran's I (spatial autocorrelation, gene 1): {moran_i:.3f}")

results = {
    "grid_size": grid_size,
    "n_spots": n_spots,
    "spatial_classification_accuracy": round(float(accuracy), 4),
    "zone_distribution": {
        "periphery": int(sum(1 for z in true_zones if z==0)),
        "intermediate": int(sum(1 for z in true_zones if z==1)),
        "core": int(sum(1 for z in true_zones if z==2))
    },
    "moran_i_gene_1": round(float(moran_i), 4),
    "interpretation": "Positive Moran\'s I indicates spatial clustering of gene expression, consistent with zonal muscle architecture",
    "recommendation": "Use spatial transcriptomics (Visium, Stereo-seq) to validate that the 30-gene panel captures zonal architecture in 3D muscle constructs"
}

(OUT / "spatial_transcriptomics.json").write_text(json.dumps(results, indent=2))
print("\nSaved to spatial_transcriptomics.json")
print("DONE")
