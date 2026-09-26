# Model card: 30-gene expression panel classifier (methods demonstration)

## Model Overview
**Name:** 30-Gene Expression Panel State Classifier  
**Version:** 1.0.0  
**Date:** 2026-05-15  
**Developer:** barlowa124 (developed independently while affiliated with the Rao Lab, NC State. The lab is not an author)  
**License:** MIT  
**Repository:** https://github.com/barlowa124/cultivated-meat-multiomic

## Intended Use
- **Primary use:** Methods demonstration. Predict cluster labels over 30 normalized gene-expression values and show the pipeline runs end to end
- **States:** expansion_competent, committed, terminal (k-means cluster labels)
- **Users:** Anyone evaluating the pipeline methodology
- **Out-of-scope:** Cultivated-meat QC, manufacturing-readiness assessment, any manufacturing decision, human diagnostics, genetic screening, food-safety pathogen detection

## Model Architecture
- **Type:** Logistic Regression (primary) + Deep Neural Network (comparison)
- **Input:** 30 normalized gene expression values (log TPM)
- **Output:** Class label + probability distribution over 3 states + confidence score
- **DNN architecture:** 30 → 64 (ReLU, dropout 0.2) → 32 (ReLU, dropout 0.2) → 3 (softmax)
- **Training:** 5-fold stratified cross-validation

## Training Data
- **Source:** Pseudo-bulk TPM + METAFlux metabolic flux derived from public human melanoma scRNA-seq datasets GSE115978 (Jerby-Arnon et al. 2018) and GSE72056 (Tirosh et al. 2016). 239 pseudo-bulk profiles used as a stand-in expression dataset
- **Labels:** K-means clustering (k=3) on PCA-reduced transcriptome + fluxome
- **Cluster sizes:** expansion_competent (61), committed (86), terminal (92)
- **Preprocessing:** log1p transformation, variance filtering (top 75%), StandardScaler

## Performance
| Metric | Value |
|--------|-------|
| CV Accuracy (Logistic) | 96.7% ± 1.0% |
| CV Accuracy (DNN) | 96.2% ± 1.5% |
| Bootstrap Accuracy | 96.3% ± 0.8% |
| Bayesian Confidence | 97.9% |
| Bovine cross-species comparison | 3-cluster structure recovered (GSE173199) |
| Noise robustness (25% CV) | 94.2% |

All accuracy figures describe classification of melanoma-derived pseudo-bulk profiles, not cultivated-meat cell states. CV figures are cross-validated; the held-out test benchmark (`Table_S2_ml_benchmark.csv`) shows the DNN dropping to 75% test accuracy where LR/RF/XGB hold ~96%. The DNN's CV number overstates its held-out performance.

## Ethical considerations & limitations
- **Domain mismatch:** Trained on human melanoma data. The panel is not a usable cultivated-meat QC assay
- **Platform limitation:** Trained on pseudo-bulk RNA-seq. No qPCR validation exists
- **Batch effects:** Models trained on a single dataset. External validation required
- **Bias:** Training data may not represent any particular cell line, media formulation, or bioreactor type

## Explainability
- SHAP top genes by rank: LMNA, PLOD1, C1QBP (per `Table_S1_panel_genes.csv` SHAP_Rank)
- Linear model coefficients provide direct interpretability

## Deployment
- **API:** Flask REST endpoint at `/predict` (see `api/app.py`)
- **Dashboard:** Streamlit interactive interface (see `dashboard/app.py`)
- **Docker:** `docker build -t cultivated-meat-multiomic .`
- **Artifacts:** `api/model.pkl`, `api/scaler.pkl`, `api/model_metadata.json`

## Maintenance
- This is a demonstration artifact. No retraining or monitoring commitments apply

## Contact
- **Issues:** https://github.com/barlowa124/cultivated-meat-multiomic/issues
