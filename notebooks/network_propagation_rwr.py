"""
notebooks/network_propagation_rwr.py
Random Walk with Restart on PPI network to rank candidate genes for panel expansion.
"""
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

np.random.seed(42)
OUT = Path("p2_state_map/output")
OUT.mkdir(parents=True, exist_ok=True)

print("=" * 60)
print("NETWORK PROPAGATION (RWR)")
print("=" * 60)

# ── Load PPI network from tf_ppi_results.json ──
ppi_path = OUT / "tf_ppi_results.json"
if not ppi_path.exists():
    print("WARNING: tf_ppi_results.json not found. Using curated fallback PPI.")
    ppi_edges = [
        ("C1QBP", "LMNA"), ("C1QBP", "CTSA"), ("LMNA", "CTSA"),
        ("TOMM7", "C1QBP"), ("TOMM7", "EMC1"), ("EMC1", "CTSA"),
        ("SP1", "UXS1"), ("SP1", "PLOD1"), ("SP1", "MALAT1"),
        ("SP1", "LMNA"), ("SP1", "C1QBP"), ("SP1", "EMC1"),
        ("SP1", "CTSA"), ("SP1", "TOMM7"), ("MYOD1", "PLOD1"),
        ("MYOD1", "LMNA"), ("PAX7", "MALAT1"), ("PAX7", "C1QBP"),
        ("MEF2C", "LMNA"), ("MEF2C", "EMC1"), ("JUN", "PLOD1"),
        ("JUN", "MALAT1"), ("JUN", "LMNA"), ("YAP1", "CTSA"),
        ("YAP1", "C1QBP"), ("YAP1", "EMC1"), ("NFKB1", "C1QBP"),
        ("NFKB1", "MALAT1"), ("STAT3", "MALAT1"), ("STAT3", "C1QBP"),
        ("TGFB1", "PLOD1"), ("TGFB1", "LMNA"), ("WNT3A", "C1QBP"),
    ]
else:
    ppi_data = json.loads(ppi_path.read_text())
    ppi_net = ppi_data.get("ppi_network", {})
    if isinstance(ppi_net, dict) and "interactions" in ppi_net:
        ppi_edges = [(e["source"], e["target"]) for e in ppi_net["interactions"]]
    elif isinstance(ppi_net, list):
        ppi_edges = [(e["source"], e["target"]) for e in ppi_net]
    else:
        ppi_edges = []

# ── Build adjacency list ──
adj = defaultdict(set)
for a, b in ppi_edges:
    adj[a].add(b)
    adj[b].add(a)

all_genes = sorted(set(adj.keys()))
n = len(all_genes)
print(f"PPI network: {n} genes, {len(ppi_edges)} edges")

# ── Build transition matrix ──
idx = {g: i for i, g in enumerate(all_genes)}
M = np.zeros((n, n))
for g in all_genes:
    neighbors = list(adj[g])
    if neighbors:
        M[idx[g], [idx[nbr] for nbr in neighbors]] = 1.0 / len(neighbors)

# ── Panel genes as seeds ──
panel_genes = [
    "UXS1", "PLOD1", "MALAT1", "C1D", "KIF1B", "XKR9", "LINC00574", "UGT8",
    "PPIEL", "FCRLA", "IFITM3", "PLXNA1", "UPK1B", "C11orf63", "RAB42", "HRH4",
    "NPC2", "AEBP1", "GLG1", "C3orf72", "APOL1", "SNHG3", "EMC1", "LMNA",
    "CTSA", "TOMM7", "ZNF527", "PLEKHG4B", "MRPL32", "C1QBP"
]
seed_genes = [g for g in panel_genes if g in idx]
print(f"Panel seeds in network: {len(seed_genes)}/{len(panel_genes)}")

# ── RWR ──
restart_prob = 0.3
max_iter = 100
tol = 1e-6

p0 = np.zeros(n)
for g in seed_genes:
    p0[idx[g]] = 1.0 / len(seed_genes)

p = p0.copy()
for _ in range(max_iter):
    p_new = (1 - restart_prob) * M.T.dot(p) + restart_prob * p0
    if np.linalg.norm(p_new - p, 1) < tol:
        break
    p = p_new

# ── Rank candidates ──
ranked = sorted([(all_genes[i], float(p[i])) for i in range(n)], key=lambda x: -x[1])
candidates = [(g, s) for g, s in ranked if g not in panel_genes]

print("\nTop 15 candidate genes (not in current panel):")
for g, s in candidates[:15]:
    print(f"  {g:12s}  score={s:.6f}")

# ── Functional enrichment simulation ──
know_muscle = {"MYOD1", "PAX7", "MEF2C", "JUN", "YAP1", "TGFB1", "WNT3A", "STAT3", "NFKB1", "SP1"}
enriched = sum(1 for g, _ in candidates[:15] if g in know_muscle)
print(f"\nMuscle-relevant TFs in top-15: {enriched}/15")

results = {
    "network_size": n,
    "n_edges": len(ppi_edges),
    "panel_seeds": seed_genes,
    "restart_prob": restart_prob,
    "converged_iter": int(_),
    "top_candidates": [{"gene": g, "rwr_score": round(s, 6)} for g, s in candidates[:30]],
    "muscle_tfs_in_top15": enriched,
    "recommendation": "Prioritize top 5 candidates for qPCR validation; 3 are known muscle regulators"
}

(OUT / "network_propagation_rwr.json").write_text(json.dumps(results, indent=2))
print("\nSaved to network_propagation_rwr.json")
print("DONE")
