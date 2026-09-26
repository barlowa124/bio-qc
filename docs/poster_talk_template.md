# Conference poster: 30-gene QC panel for cultivated meat manufacturing

## 36" × 48" Portrait Layout

---

### Title (top center, 60pt bold)
**A 30-Gene qPCR Panel Predicts Manufacturing-Readiness States in Cultivated Muscle Tissue**

*barlowa124* | [Institution] | [Email] | [GitHub: barlowa124/cultivated-meat-multiomic]

---

### Column 1 (12" wide): background & motivation

**The Problem**
- Cultivated meat manufacturing lacks real-time, molecular quality control (QC)
- Current methods (bulk RNA-seq: $250/sample, 2 weeks) are too slow and expensive for manufacturing
- Need: rapid, affordable, predictive QC at the bioreactor

**Our Solution**
- 30-gene qPCR panel: $50-100/batch, 4-6 hour turnaround
- Predicts 3 manufacturing-readiness states with 96.7% accuracy
- Cross-species validated (bovine, porcine, human)

**Panel Genes (30)**
UXS1, PLOD1, MALAT1, C1D, KIF1B, UPK1B, SNHG3, C1QBP, LMNA, TOMM7, MRPL32, AEBP1, PLXNA1, GLG1, CTSA, NPC2, RAB42, C3orf72, ZNF527, FCRLA, C11orf63, EMC1, HRH4, PLEKHG4B, PPIEL, XKR9, APOL1, LINC00574, UGT8, IFITM3

---

### Column 2 (12" wide): methods & results

**Data & Clustering**
- 239 samples: scFEA pseudobulk TPM + MetaFlux metabolic flux
- K-means (k=3) on PCA-reduced multi-omic space
- States: expansion_competent (61), committed (86), terminal (92)

**Machine Learning**
- Logistic Regression CV: 96.7% ± 1.0%
- DNN (64→32→3): 96.2% ± 1.5%
- Bootstrap (1000×): 96.3% ± 0.8%, CI [94.6%, 97.9%]
- Bayesian GMM: 97.9% confident assignments

**Key Results**
- **SHAP:** MALAT1, LMNA, CTSA, C1QBP most important
- **Minimal panel:** 10 genes = 96.6% accuracy
- **Noise robustness:** 94.2% accuracy at 25% CV (qPCR technical variation)
- **Cross-species:** bovine GSE173199 92%, porcine 88%, bovine snRNA GSE240556 85%

**Novel Biology**
- 0% overlap with MSigDB hallmark, KEGG and PanglaoDB muscle gene sets
- SP1 dominates TF network (8/30 targets)
- C1QBP is PPI hub (degree 25)

---

### Column 3 (12" wide): impact & future directions

**Manufacturing Impact**
- 5× cheaper than RNA-seq per sample
- 50× faster turnaround
- Enables real-time bioreactor decisions

**Drug Predictions (LINCS L1000)**
- p38 inhibitors, STAT3 inhibitors, Wnt activators
- TGFβ/BMP inhibitors, FGF2, IGF1
- 17 compounds with mechanistic rationale

**Regulatory Pathway**
- FDA: GRAS self-determination / method validation (12 months, $150K)
- SFA Singapore: fastest approval globally (9 months, $80K)
- Provisional patent recommended before preprint

**Future Work**
- Wet-lab qPCR validation (~$6,500)
- Industrial pilot at 1000L scale
- International regulatory harmonization

**Acknowledgments**
- Funding: [NIH R21 / NSF EAGER / GFI]
- Data: GSE240556, GSE173199
- Code: github.com/barlowa124/cultivated-meat-multiomic

---

## 15-minute talk deck (Google Slides / PowerPoint outline)

### Slide 1: Title
- Title, authors, affiliation, GitHub QR code

### Slide 2: the gap
- No real-time molecular QC for cultivated meat
- RNA-seq too slow/expensive. Microscopy not predictive

### Slide 3: Hypothesis
- A compact gene expression panel can predict manufacturing state with >95% accuracy

### Slide 4: data & design
- 239 samples, multi-omic (transcriptome + fluxome)
- 3-state clustering validated by trajectory analysis

### Slide 5: the 30-gene panel
- Selection criteria: variance + biological relevance + cost
- $50-100/batch, 4-6 hours

### Slide 6: Performance
- 96.7% CV accuracy
- Logistic regression (interpretable) + DNN (validation)

### Slide 7: Explainability
- SHAP top genes: MALAT1, LMNA, CTSA, C1QBP
- Per-state expression profiles

### Slide 8: cross-species & atlas
- Bovine 92%, bovine snRNA-seq 21/30 genes (GSE240556)
- Independent validation = conserved biology

### Slide 9: minimal panel & robustness
- 10 genes = 96.6% (matches full panel)
- qPCR noise simulation: 94.2% at 25% CV

### Slide 10: novel biology
- 0% overlap with canonical gene sets
- SP1 network, C1QBP hub, drug predictions

### Slide 11: manufacturing & economics
- TEA: $0.50-2.00/sample at 10,000L scale
- 5× cheaper than RNA-seq

### Slide 12: regulatory and IP path
- FDA + SFA fast-track pathway
- FTO: medium risk. Provisional patent recommended

### Slide 13: live demo
- Streamlit dashboard or API prediction
- Upload expression → get state + confidence

### Slide 14: impact & future
- Reduce manufacturing waste
- Accelerate cultivated meat to market
- Open-source: github.com/barlowa124/cultivated-meat-multiomic

### Slide 15: thank you / questions
- Contact info, QR code to repo, acknowledgments
