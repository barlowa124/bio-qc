# Multi-omic state-map pipeline: RNA + metabolic-flux clustering with cross-species validation

[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

**Methods demonstration on public datasets. Developed independently by the author while affiliated with the Rao Lab, NC State. See Scope.**

---

## Scope and data provenance (read first)

The state-map and gene-panel analyses in `p2_state_map/` and `p3_qc_panel/` were developed on pseudo-bulk profiles derived from two public human melanoma single-cell RNA-seq datasets (GSE115978, Jerby-Arnon et al. 2018, and GSE72056, Tirosh et al. 2016). These were used as a stand-in expression dataset while building the pipeline. The resulting clusters, accuracy figures and gene rankings therefore characterize melanoma-derived profiles and are **not** evidence about cultivated-meat cell states, manufacturing readiness, or a usable qPCR QC panel.

The bovine (GSE173199), porcine (GSE206914) and bovine single-nucleus (GSE240556) analyses use skeletal-muscle/myogenic data and are the domain-relevant components of this repository.

Downstream artifacts generated from the melanoma-derived panel (drug-response scoring, techno-economic, life-cycle, regulatory and patent drafts) are illustrative pipeline outputs demonstrating that the tooling runs end to end. They are not findings and should not be cited as such.

This repository was originally framed as a cultivated meat manufacturing-QC project. It is now maintained as a methods demonstration.

## Pipeline demonstration metrics

| Metric | Value |
|--------|-------|
| State-map cluster separation accuracy (melanoma pseudo-bulk) | **96.7%** (5-fold CV) |
| 30-gene panel accuracy (melanoma pseudo-bulk) | **96.7%** ± 1.0% |
| Bootstrap stability, 1000 iterations (melanoma pseudo-bulk) | **96.3%** ± 0.8%, 95% CI [94.6%, 97.9%] |
| Bayesian confident assignments (melanoma pseudo-bulk) | **97.9%** |
| Cross-species comparison | 3-cluster structure recovered in bovine GSE173199 |
| snRNA-seq comparison | 17,541 bovine muscle nuclei (GSE240556), 21/30 panel genes detected |
| PPI network | 74 edges, C1QBP central hub (degree 25) |
| Top TF regulator | SP1 (8/30 targets, 26.7%) |

Accuracy figures measure how well a classifier recovers k-means cluster labels that were themselves defined on the same expression embedding. They quantify cluster separability and panel sufficiency, not agreement with any external ground truth or biological state annotation.

## Three expression clusters (k-means, k=3)

The pipeline partitions the melanoma-derived pseudo-bulk embedding into three clusters, labeled by their marker profiles:

1. **expansion_competent**: high proliferative-capacity markers, low differentiation markers
2. **committed**: intermediate metabolic activity, mixed marker expression
3. **terminal**: high differentiation markers, low proliferative-capacity markers

## Quick Start

```bash
# Clone
git clone https://github.com/barlowa124/cultivated-meat-multiomic.git
cd cultivated-meat-multiomic

# Install
pip install -r requirements.txt

# Run core analysis (requires the source expression/flux matrices —
# see docs/installation.md; analysis scripts read shared_data/data_paths.yaml
# or absolute paths from the original analysis workstation)
python notebooks/comprehensive_analysis.py

# Run high-value analyses
python notebooks/porcine_3species.py        # Cross-species comparison
python notebooks/snrna_seq_analysis.py      # Single-nucleus comparison
python notebooks/vae_bayesian.py            # VAE + Bayesian GMM
python notebooks/tf_ppi_analysis.py         # TF enrichment + PPI network
python notebooks/bootstrap_ml_timeseries.py # Bootstrap + ML + time-series
python notebooks/shap_dnn.py                # SHAP explainability + DNN
python notebooks/wgcna_de_batch.py          # WGCNA + DE + batch effects
python notebooks/cellcom_benchmark.py       # Cell communication + benchmark

# Generate report
python notebooks/build_report.py            # HTML report → Print to PDF

# Start prediction API
cd api && pip install -r requirements.txt && python app.py
# Binds 127.0.0.1 by default; set HOST to override.
# The API loads pickled model artifacts from this repo only; it is for
# local/research use and is not a hardened service.
```

### Docker

```bash
docker build -t cultivated-meat .
docker run -v ./data:/app/data -v ./output:/app/output cultivated-meat
```

## Project Structure

```
cultivated_meat_projects/
├── notebooks/                  # Analysis scripts
│   ├── comprehensive_analysis.py   # Master pipeline
│   ├── p2_state_map.py             # Multi-omic state map
│   ├── p3_qc_panel.py              # 30-gene panel selection
│   ├── shap_dnn.py                 # SHAP + deep learning
│   ├── wgcna_de_batch.py           # Co-expression + DE
│   ├── cellcom_benchmark.py        # Cell communication + benchmark
│   └── ...                         # many more scripts
├── api/                       # Flask prediction API
│   ├── app.py                     # Server
│   ├── model.pkl                  # Trained model
│   └── test_client.py             # Test client
├── docs/                      # Jupyter Book documentation
├── cross_species_validation/  # Bovine/porcine GEO data
├── p2_state_map/output/       # State-map results, figures, reports
├── p3_qc_panel/output/        # 30-gene panel results
├── Dockerfile                 # Container definition
└── docker-compose.yml         # Multi-service orchestration
```

## Analyses Included

| Analysis | Description |
|----------|-------------|
| P2 State Map | Joint RNA+flux PCA + K-means clustering |
| P3 Panel | L1 logistic regression biomarker selection |
| Cross-Species | Bovine (GSE173199) + Porcine (GSE206914) comparison |
| snRNA-seq | GSE240556 single-nucleus resolution (17,541 nuclei) |
| VAE + Bayesian GMM | Deep embeddings + uncertainty quantification |
| TF Enrichment | JASPAR-based transcription factor analysis (30 TFs) |
| PPI Network | STRING protein-protein interaction network (74 edges) |
| SHAP Explainability | Per-sample feature importance |
| Bootstrap Stability | 1000-iteration gene selection stability |
| ML Comparison | LR vs RF vs SVM vs DNN benchmark |
| Time-Series | Bovine D0-D7 temporal gene trends |
| Cell Communication | Ligand-receptor analysis from snRNA-seq |
| Integration Benchmark | Early vs late vs single-modality comparison |
| Differential Expression | 75 significant state comparisons |
| Batch Effects | Cross-batch variability assessment |
| Batch Correction | ComBat/Harmony/MNN/reference atlas benchmarking |
| Pathway Enrichment | GO BP / KEGG / Reactome (hypergeometric, BH-corrected) |
| Cross-Platform | RNA-seq vs qPCR vs Nanostring concordance |
| Transfer Learning | Human→porcine domain adaptation (+8.3% accuracy) |
| Conformal Prediction | 91.2% coverage, 78% singleton sets |
| Adversarial Robustness | FGSM perturbation (epsilon 0.01–0.1) stability |
| Pareto Optimization | Cost-accuracy-gene count tradeoff front |
| CRISPR Screen Mining | Panel-overlap workflow on *simulated* DepMap-style scores (`crispr_screen_mining.py` seeds np.random — wiring real DepMap/Avana data is a TODO; treat scores as demo only) |
| Alternative Splicing | Isoform ratio state discrimination (12/30 genes) |
| Reference Atlas | Human Protein Atlas outlier flagging |
| Concept Drift | KS-test + Mahalanobis monthly retraining triggers |

Additional exploratory scripts (illustrative only): `patent_claims.py`, `patent_landscape.py`, `regulatory_dossier.py`, `regulatory_pathway.py`, `techno_economic_analysis.py`, `commercial_qc_benchmark.py`, `supply_chain_risk.py`, `lca_comparison.py`, `microbiome_primers.py`, `digital_twin.py`, `multiomics_integration.py` cost comparisons. These draft techno-economic, regulatory, patent, supply-chain, life-cycle, microbiome-screen, digital-twin and cost outputs as end-to-end pipeline demonstrations. They are not findings.

## Interactive figures

| Figure | Content | File |
|--------|---------|------|
| Fig 1 | State Map (UMAP scatter) | `docs/figures/fig1_state_map.html` |
| Fig 2 | Panel ML Benchmark (6 classifiers) | `docs/figures/fig2_qc_performance.html` |
| Fig 3 | SHAP Gene Importance (top 15) | `docs/figures/fig3_shap_importance.html` |
| Fig 4 | Cross-Species Comparison (accuracy + sample counts) | `docs/figures/fig4_cross_species.html` |
| Fig 5 | Noise Robustness (CV vs accuracy) | `docs/figures/fig5_noise_robustness.html` |
| Fig 6 | Drug Prediction Scores (LINCS/CMap) | `docs/figures/fig6_drug_predictions.html` |
| Fig 7 | Techno-Economic Sensitivity | `docs/figures/fig7_tea_sensitivity.html` |
| Fig 9 | Batch Correction Benchmark | `docs/figures/fig9_batch_correction.html` |
| Fig 10 | Pathway Enrichment Dot Plot | `docs/figures/fig10_pathway_enrichment.html` |
| Fig 11 | Cross-Platform Comparison | `docs/figures/fig11_cross_platform.html` |
| Fig 11b | Platform Heatmap (20 genes x 3 platforms) | `docs/figures/fig11b_platform_heatmap.html` |

All figures are standalone HTML. Double-click to open or embed in presentations. Each file loads Plotly.js v3.3.1 from the CDN (`cdn.plot.ly`), so an internet connection is required for interactive rendering.

**Generate:** `python notebooks/generate_interactive_figures.py`

## ML rigor analyses

| Analysis | Key Result | File |
|----------|-----------|------|
| Ensemble (LR+DNN+XGB) | **97.9%** test accuracy | `docs/figures/figA1_ensemble.svg` |
| Calibration Curves | Brier scores + reliability diagrams | `docs/figures/figA2_calibration.svg` |
| Multi-class ROC/PR | ROC-AUC = **0.9967**, AP = **0.9950** | `docs/figures/figA3_roc_pr.svg` |
| Ablation Study | MALAT1 drop **12.5%**, LMNA drop **10.4%** | `docs/figures/figA4_ablation.svg` |
| Learning Curves | Plateau by ~120 samples | `docs/figures/figA5_learning_curves.svg` |
| Bootstrap SHAP Stability | LMNA, TOMM7, C1QBP most stable (CV < 0.11) | `docs/figures/figA6_bootstrap_shap.svg` |

**Generate:** `python notebooks/ml_rigor_analyses.py`

## Supplementary Tables

All supplementary tables are generated programmatically from JSON outputs:

| Table | Content | File |
|-------|---------|------|
| S1 | Full 30-gene panel (gene, function, SHAP rank, primer) | `docs/supplementary/Table_S1_panel_genes.csv` |
| S2 | ML benchmark summary (6 classifiers + ensemble) | `docs/supplementary/Table_S2_ml_benchmark.csv` |
| S3 | Cross-species concordance (bovine, porcine, human) | `docs/supplementary/Table_S3_cross_species.csv` |
| S4 | Cross-platform correlation matrix | `docs/supplementary/Table_S4_cross_platform.csv` |
| S5 | Batch correction comparison (5 methods) | `docs/supplementary/Table_S5_batch_correction.csv` |
| S6 | Pathway enrichment results (GO/KEGG/Reactome) | `docs/supplementary/Table_S6_pathway_enrichment.csv` |
| S7 | Bootstrap SHAP stability rankings (100 resamples) | `docs/supplementary/Table_S7_bootstrap_shap.csv` |
| S8 | Ablation study (leave-one-gene-out drops) | `docs/supplementary/Table_S8_ablation.csv` |

**Generate:** `python notebooks/generate_supplementary_tables.py`

## Data Sources

- **RNA-seq**: 239 pseudo-bulk profiles derived from human melanoma scRNA-seq (GSE115978 and GSE72056), 23,682 genes, scFEA/METAFlux flux estimates
- **METAFlux**: 13,082 reactions aggregated to 10 metabolic pathways
- **Bovine**: GSE173199 (38 samples, D0-D7 timecourse)
- **Porcine**: GSE206914 (45 samples, embryonic stages)
- **snRNA-seq**: GSE240556 (17,541 bovine muscle nuclei)

## Acknowledgements

Developed independently by barlowa124 while affiliated with the Rao Lab, North Carolina State University. Public data from GEO as listed above. The lab is not an author of this repository.

## License

MIT License. See [LICENSE](LICENSE) file for details.
