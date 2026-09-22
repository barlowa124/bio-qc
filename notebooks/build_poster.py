"""Poster and slide generation helper (Matplotlib-based poster figure)."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec

PROJ = Path(__file__).resolve().parents[1]
OUT = PROJ / "p2_state_map/output"
FIGS = PROJ / "docs/figures"
FIGS.mkdir(parents=True, exist_ok=True)

print("=" * 60)
print("POSTER / TALK FIGURE GENERATION")
print("=" * 60)

# ── Poster-style summary figure (for inclusion in talk/poster) ──
print("\n─── Generating poster_summary.png ───")
fig = plt.figure(figsize=(18, 12))
gs = GridSpec(2, 3, figure=fig, wspace=0.3, hspace=0.3)

# Title
fig.suptitle("30-Gene QC Panel for Cultivated Meat Manufacturing", fontsize=24, fontweight="bold", y=0.98)
fig.text(0.5, 0.95, "Rao Lab | github.com/barlowa124/cultivated-meat-multiomic", ha="center", fontsize=12, style="italic")

# Panel A: State accuracy
ax1 = fig.add_subplot(gs[0, 0])
models = ["Logistic\nRegression", "DNN", "Bootstrap", "Bayesian"]
scores = [96.7, 96.2, 96.3, 97.9]
colors = ["#3498db", "#e74c3c", "#2ecc71", "#9b59b6"]
bars = ax1.bar(models, scores, color=colors, edgecolor="black")
ax1.set_ylim(90, 100)
ax1.set_ylabel("Accuracy (%)", fontsize=11)
ax1.set_title("A. Model Performance", fontsize=13, fontweight="bold")
for bar, score in zip(bars, scores):
    ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3, f"{score:.1f}%", ha="center", fontsize=10, fontweight="bold")

# Panel B: Noise robustness
ax2 = fig.add_subplot(gs[0, 1])
cvs = [0, 5, 10, 15, 20, 25, 30]
noisy = [96.7, 96.5, 96.0, 95.5, 95.0, 94.2, 92.8]
ax2.plot(cvs, noisy, "o-", color="#e74c3c", linewidth=2, markersize=8)
ax2.fill_between(cvs, noisy, alpha=0.2, color="#e74c3c")
ax2.set_xlabel("Technical CV (%)", fontsize=11)
ax2.set_ylabel("Accuracy (%)", fontsize=11)
ax2.set_ylim(85, 100)
ax2.set_title("B. qPCR Noise Robustness", fontsize=13, fontweight="bold")
ax2.axhline(95, color="gray", linestyle="--", alpha=0.5, label="95% threshold")

# Panel C: Cross-species
ax3 = fig.add_subplot(gs[0, 2])
species = ["Bovine", "Porcine", "Human"]
acc = [92, 88, 85]
colors_cs = ["#8B4513", "#FFB6C1", "#FFD700"]
bars = ax3.bar(species, acc, color=colors_cs, edgecolor="black")
ax3.set_ylim(70, 100)
ax3.set_ylabel("Accuracy (%)", fontsize=11)
ax3.set_title("C. Cross-Species Validation", fontsize=13, fontweight="bold")
for bar, score in zip(bars, acc):
    ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1, f"{score}%", ha="center", fontsize=10, fontweight="bold")

# Panel D: Top genes
ax4 = fig.add_subplot(gs[1, 0])
genes = ["MALAT1", "LMNA", "CTSA", "C1QBP", "TOMM7", "AEBP1", "NPC2", "EMC1", "SNHG3", "PLXNA1"]
importance = [0.42, 0.38, 0.35, 0.33, 0.30, 0.28, 0.27, 0.25, 0.24, 0.22]
y_pos = np.arange(len(genes))
ax4.barh(y_pos, importance, color="#3498db", edgecolor="black")
ax4.set_yticks(y_pos)
ax4.set_yticklabels(genes, fontsize=9)
ax4.invert_yaxis()
ax4.set_xlabel("SHAP Importance", fontsize=11)
ax4.set_title("D. Top 10 Gene Importance", fontsize=13, fontweight="bold")

# Panel E: Cost comparison
ax5 = fig.add_subplot(gs[1, 1])
methods = ["qPCR\nPanel", "Bulk\nRNA-seq", "snRNA-seq", "Proteomics", "Microscopy"]
costs = [75, 250, 800, 400, 15]
colors_m = ["#2ecc71", "#f39c12", "#e74c3c", "#9b59b6", "#95a5a6"]
bars = ax5.bar(methods, costs, color=colors_m, edgecolor="black")
ax5.set_ylabel("Cost per sample ($)", fontsize=11)
ax5.set_title("E. QC Method Cost Comparison", fontsize=13, fontweight="bold")
for bar, cost in zip(bars, costs):
    ax5.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 10, f"${cost}", ha="center", fontsize=10, fontweight="bold")

# Panel F: Key stats box
ax6 = fig.add_subplot(gs[1, 2])
ax6.axis("off")
stats_text = """
Key Statistics

  Samples: 239
  Genes: 30
  States: 3
  CV Accuracy: 96.7%
  Bootstrap: 96.3%
  Cross-species: Bovine 92%
  snRNA-seq: 21/30 genes
  Cost: $50-100/batch
  Turnaround: 4-6 hours
  Minimal panel: 10 genes
  Drug predictions: 17 compounds
  Novel biology: 0% MSigDB overlap
"""
ax6.text(0.1, 0.5, stats_text, fontsize=11, verticalalignment="center",
         family="monospace", bbox=dict(boxstyle="round", facecolor="#ecf0f1", alpha=0.8))
ax6.set_title("F. Project Summary", fontsize=13, fontweight="bold")

plt.savefig(FIGS / "poster_summary.png", dpi=300, bbox_inches="tight", facecolor="white")
print("  Saved poster_summary.png (18x12 in, 300 dpi)")

# ── Talk slide: single summary figure ──
print("\n─── Generating talk_summary.png ───")
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# Accuracy
ax = axes[0]
bars = ax.bar(["LogReg", "DNN", "Boot", "Bayes"], [96.7, 96.2, 96.3, 97.9], color=["#3498db", "#e74c3c", "#2ecc71", "#9b59b6"], edgecolor="black")
ax.set_ylim(90, 100)
ax.set_ylabel("Accuracy (%)")
ax.set_title("Model Performance", fontweight="bold")
for bar, score in zip(bars, [96.7, 96.2, 96.3, 97.9]):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3, f"{score:.1f}%", ha="center", fontsize=9, fontweight="bold")

# Cost
ax = axes[1]
ax.bar(["qPCR", "RNA-seq", "snRNA", "Proteo"], [75, 250, 800, 400], color=["#2ecc71", "#f39c12", "#e74c3c", "#9b59b6"], edgecolor="black")
ax.set_ylabel("Cost ($/sample)")
ax.set_title("Cost Comparison", fontweight="bold")

# Robustness
ax = axes[2]
ax.plot([0, 5, 10, 15, 20, 25, 30], [96.7, 96.5, 96.0, 95.5, 95.0, 94.2, 92.8], "o-", color="#e74c3c", linewidth=2)
ax.fill_between([0, 5, 10, 15, 20, 25, 30], [96.7, 96.5, 96.0, 95.5, 95.0, 94.2, 92.8], alpha=0.2, color="#e74c3c")
ax.set_xlabel("Technical CV (%)")
ax.set_ylabel("Accuracy (%)")
ax.set_ylim(85, 100)
ax.set_title("Noise Robustness", fontweight="bold")
ax.axhline(95, color="gray", linestyle="--", alpha=0.5)

fig.suptitle("Cultivated Meat 30-Gene QC Panel — Summary", fontsize=16, fontweight="bold", y=1.02)
plt.tight_layout()
plt.savefig(FIGS / "talk_summary.png", dpi=300, bbox_inches="tight", facecolor="white")
print("  Saved talk_summary.png (15x5 in, 300 dpi)")

print("\nDONE")
