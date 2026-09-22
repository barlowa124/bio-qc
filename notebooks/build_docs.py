"""Build Jupyter Book documentation site."""
from pathlib import Path

PROJ = Path(__file__).resolve().parents[1]
DOCS = PROJ / "docs"
DOCS.mkdir(exist_ok=True)

print("=" * 60)
print("JUPYTER BOOK DOCUMENTATION")
print("=" * 60)

# _config.yml
config = """# Jupyter Book configuration for Cultivated Meat Analysis Pipeline
title: "Cultivated Meat Multi-Omic Analysis"
author: "Rao Lab, NC State University"
logo: ""
execute:
  execute_notebooks: "off"
repository:
  url: https://github.com/rao-lab/cultivated-meat-multiomic
  branch: main
html:
  use_issues_button: true
  use_repository_button: true
  home_page_in_navbar: true
launch_buttons:
  colab_url: "https://colab.research.google.com"
parse:
  myst_enable_extensions:
    - dollarmath
    - linkify
    - substitution
sphinx:
  config:
    html_theme: sphinx_book_theme
    html_theme_options:
      show_toc_level: 2
"""

(DOCS / "_config.yml").write_text(config)

# _toc.yml
toc = """format: jb-book
root: index
parts:
  - caption: Overview
    chapters:
    - file: overview
    - file: installation
  - caption: Core Analyses
    chapters:
    - file: p2_state_map
    - file: p3_qc_panel
    - file: p5_reproducibility
  - caption: Extended Analyses
    chapters:
    - file: p1_protocol_meta
    - file: p4_media_formulation
    - file: cross_species
    - file: snrna_seq
    - file: vae_bayesian
    - file: tf_ppi
    - file: bootstrap_ml
  - caption: Results
    chapters:
    - file: figures
    - file: manuscript
  - caption: Appendix
    chapters:
    - file: supplementary
    - file: validation_plan
"""

(DOCS / "_toc.yml").write_text(toc)

# index.md
index_md = """# Cultivated Meat Multi-Omic Analysis Pipeline

## A Manufacturing-Readiness State Map with 30-Gene qPCR QC Panel

**Rao Lab, North Carolina State University**

---

### Key Results

| Metric | Value |
|--------|-------|
| Multi-omic state map accuracy | **96.7%** (5-fold CV) |
| 30-gene panel accuracy | **96.7%** ± 1.0% |
| Cross-species validation | Bovine 3-state conservation confirmed |
| snRNA-seq validation | 17,541 nuclei, 21/30 panel genes detected |
| Bayesian confidence | 97.9% confident assignments |
| Assay cost | \\$50-100/batch |
| Turnaround time | 4-6 hours |

### Quick Start

```bash
# Clone
git clone https://github.com/rao-lab/cultivated-meat-multiomic
cd cultivated-meat-multiomic

# Install
pip install -r requirements.txt

# Run core analysis
python notebooks/comprehensive_analysis.py

# Run high-value analyses
python notebooks/porcine_3species.py
python notebooks/snrna_seq_analysis.py
python notebooks/vae_bayesian.py
python notebooks/tf_ppi_analysis.py
python notebooks/bootstrap_ml_timeseries.py
```

### Docker

```bash
docker build -t cultivated-meat .
docker run -v ./data:/app/data -v ./output:/app/output cultivated-meat
```

### Citation

If you use this work, please cite our preprint (forthcoming on bioRxiv).
"""

(DOCS / "index.md").write_text(index_md)

# overview.md
overview_md = """# Project Overview

## Motivation

Cultivated meat production requires robust, quantitative quality control metrics to assess stem cell manufacturing readiness. Current approaches rely on subjective morphological assessment and endpoint assays that cannot provide real-time feedback for bioprocess optimization.

## Approach

We integrate:
- **RNA-seq transcriptomics** (425 samples across conditions)
- **METAFlux metabolic flux predictions** (239 samples with paired flux)
- **Cross-species validation** (bovine GSE173199, GSE240556; porcine GSE206914)
- **Machine learning** (PCA, K-means, L1 logistic regression, VAE, Bayesian GMM)

## Three Manufacturing Readiness States

1. **Expansion-competent** — High proliferative capacity, low differentiation markers
2. **Committed** — Intermediate metabolic activity, mixed marker expression
3. **Terminal** — High differentiation markers, low proliferative capacity

## 30-Gene qPCR Panel

The panel achieves 96.7% cross-validation accuracy matching the full multi-omic embedding, enabling routine batch triage at \\$50-100 per assay with 4-6 hour turnaround.
"""

(DOCS / "overview.md").write_text(overview_md)

# installation.md
install_md = """# Installation

## Requirements

- Python 3.10+
- 16GB+ RAM recommended
- CUDA-capable GPU optional (for VAE training)

## Dependencies

```bash
pip install -r requirements.txt
```

Key packages:
- `numpy`, `pandas`, `scipy` — Data processing
- `scikit-learn` — Machine learning
- `torch` — Deep learning (VAE)
- `anndata`, `scipy.sparse` — Single-cell data
- `matplotlib`, `seaborn` — Visualization

## Data

Public datasets are downloaded automatically or can be placed in:
- `cross_species_validation/` — Bovine and porcine GEO data
- Source data paths configured in `shared_data/data_paths.yaml`
"""

(DOCS / "installation.md").write_text(install_md)

# Create placeholder pages for each analysis
analyses = {
    "p2_state_map": "# Multi-Omic State Map\n\nJoint PCA embedding of RNA-seq and METAFlux data with K-means clustering.",
    "p3_qc_panel": "# 30-Gene QC Panel\n\nL1-regularized logistic regression for minimal biomarker selection.",
    "p5_reproducibility": "# Reproducibility Analysis\n\nGene stability, PCA variance, and clustering stability across operators.",
    "p1_protocol_meta": "# Protocol Meta-Analysis\n\nParsing protocol documents and lab spreadsheets for factor extraction.",
    "p4_media_formulation": "# Media Formulation Analysis\n\nComponent-response signatures and optimization recommendations.",
    "cross_species": "# Cross-Species Validation\n\nBovine and porcine data processing with 3-species comparison.",
    "snrna_seq": "# snRNA-seq Analysis\n\nGSE240556 single-nucleus resolution validation of the 30-gene panel.",
    "vae_bayesian": "# VAE + Bayesian GMM\n\nVariational autoencoder embeddings and Bayesian uncertainty quantification.",
    "tf_ppi": "# TF Enrichment + PPI Network\n\nTranscription factor binding site analysis and protein-protein interaction network.",
    "bootstrap_ml": "# Bootstrap + ML Comparison\n\nGene panel stability via bootstrap and multi-model benchmarking.",
    "figures": "# Publication Figures\n\nFive publication-quality figures with consistent styling.",
    "manuscript": "# Manuscript\n\nFull manuscript with cross-species results.",
    "supplementary": "# Supplementary Materials\n\nTables S1-S4, Figures S1-S3.",
    "validation_plan": "# Validation Plan\n\nProspective lab validation protocol for 30-batch study.",
}

for name, content in analyses.items():
    (DOCS / f"{name}.md").write_text(content)

print(f"Created {len(analyses)} documentation pages")
print(f"\nDocumentation ready at: {DOCS}")
print("To build: pip install jupyter-book && jb build docs/")
print("=" * 60)
print("DONE")
