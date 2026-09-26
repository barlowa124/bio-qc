# Foundation model integration plan
## Cultivated meat 30-gene qPCR panel: 3-state classification

---

## 1. Model details

### Geneformer (ctheodoris/Geneformer)
| Property | Value |
|----------|-------|
| Architecture | BERT (encoder-only transformer) |
| Parameters | 316M (V2), also 10M (V1), 104M (V2) |
| Pretraining Data | Genecorpus-104M: ~104M human single-cell transcriptomes |
| Input Encoding | Rank-value encoding (genes ranked by expression within cell, scaled across corpus) |
| Objective | Masked language modeling (15% masking) |
| License | Apache 2.0 |
| Downloads | 4.1M |
| Key Capability | Zero-shot in silico perturbation + fine-tuning for classification |

**How it works:**
1. Each cell's transcriptome → rank genes by expression → tokenize as rank values
2. 15% of genes masked → model predicts masked gene from context
3. Attention weights learn gene network hierarchy in self-supervised manner
4. Fine-tune by replacing head with task-specific classifier

**Relevance to cultivated meat:**
- Cell state classification across differentiation. Directly matches your 3-state problem
- Batch integration. Addresses your batch correction needs
- Disease classification. Analogous to manufacturing-readiness state detection
- In silico perturbation. Could identify which genes drive state transitions

**Caveats:**
- Pretrained on human only. Cross-species transfer to bovine/porcine is unknown
- Input size 4096 genes. Your 30-gene panel is much smaller and needs gene set intersection

### scGPT (tdc/scGPT)
| Property | Value |
|----------|-------|
| Architecture | Custom generative transformer (scGPT) |
| Parameters | 50.8M |
| Pretraining Data | 33M cells from CELLxGENE |
| Input Encoding | Gene token embeddings + expression values |
| Objective | Generative pretraining (auto-regressive) |
| License | MIT |
| Downloads | 12.4K |
| Key Capability | Multi-omic integration + perturbation prediction |

**How it works:**
1. Genes treated as tokens, expression as values
2. Auto-regressive generation of gene expression
3. Specialized attention for gene-gene interactions
4. Fine-tune via transfer learning on downstream tasks

**Relevance to cultivated meat:**
- Multi-batch integration. Directly applicable
- Cell type annotation. Matches your classification task
- Multi-omic integration. Could incorporate proteomics/metabolomics later
- MIT license. Most permissive for commercial use

**Caveats:**
- Smaller model (50M). May not capture as much biology as Geneformer
- Less community adoption (12K vs 4.1M downloads)

### sCellTransformer (InstaDeepAI/sCellTransformer)
| Property | Value |
|----------|-------|
| Architecture | Custom sCT |
| Parameters | 79.2M |
| Downloads | 375 |
| License | Not specified |

**Verdict:** Too niche, low adoption. Skip for now.

---

## 2. Can they be fine-tuned for 3-state classification?

### Geneformer: yes
The README explicitly lists "cell state classification across differentiation" as a demonstrated application. The fine-tuning workflow:
1. Tokenize your 30-gene qPCR data using Geneformer's tokenizer
2. Map your genes to Geneformer's vocabulary (~20K protein-coding genes)
3. Add a classification head (3 output neurons for Proliferative/Differentiating/Mature)
4. Fine-tune with your 239 labeled samples
5. Evaluate against your existing LR/DNN/XGB benchmarks

**Challenge:** Geneformer expects 4096 genes as input. Your 30-gene panel is much smaller. Options:
- **A) Gene set intersection:** Only use the 30 genes, pad the rest with zeros. Geneformer's attention will focus on your genes.
- **B) Imputation:** Use Geneformer to impute the missing ~4000 genes from your 30, then classify on full profile.
- **C) Embedding extraction:** Extract cell embeddings from Geneformer, train a lightweight classifier on top.

### scGPT: yes
The README lists "cell type annotation" as a downstream task. The fine-tuning workflow:
1. Load your data as AnnData
2. Use TDC's scGPT tokenizer
3. Fine-tune the pretrained model
4. Simpler API than Geneformer

**Challenge:** scGPT expects full transcriptome. Same gene set intersection issue.

---

## 3. Integration Plan

### Phase 1: embedding extraction (low risk, 1-2 hours)
**Goal:** Test if foundation model embeddings improve over your current features.

1. Install Geneformer: `pip install git+https://huggingface.co/ctheodoris/Geneformer`
2. Map your 30 genes to Geneformer vocabulary
3. Extract cell embeddings for all 239 samples
4. Train LR/SVM on embeddings vs raw qPCR Ct values
5. Compare accuracy

**Success metric:** Embedding-based classifier > 0.967 (current LR baseline)

### Phase 2: full fine-tuning (medium risk, 3-5 hours)
**Goal:** Fine-tune Geneformer end-to-end on your 3-state task.

1. Prepare tokenized dataset from qPCR data
2. Add 3-class classification head
3. Fine-tune with low learning rate (1e-5)
4. Compare against your ensemble (0.979)
5. Run SHAP on fine-tuned model for gene importance

**Success metric:** Fine-tuned Geneformer > 0.979 (current ensemble)

### Phase 3: cross-species transfer (high impact, 3-5 hours)
**Goal:** Use Geneformer's human knowledge to boost bovine/porcine accuracy.

1. Fine-tune Geneformer on human data (Phase 2)
2. Map bovine/porcine orthologs to human genes
3. Evaluate zero-shot on bovine (GSE173199) and porcine (GSE206914)
4. Fine-tune further on bovine/porcine if needed

**Success metric:** Bovine accuracy > 0.92, Porcine > 0.89 (current baselines)

### Phase 4: in silico perturbation (novel contribution, 4-6 hours)
**Goal:** Use Geneformer's zero-shot perturbation to identify genes that drive state transitions.

1. For each of the 30 genes, simulate knockout/overexpression
2. Predict resulting cell state shift
3. Rank genes by perturbation impact
4. Compare with SHAP importance and ablation study rankings
5. This becomes a novel figure/table for the manuscript

**Success metric:** Perturbation rankings correlate with SHAP (Spearman ρ > 0.5)

---

## 4. Risk Assessment

| Risk | Likelihood | Impact | Mitigation |
|------|-----------|--------|------------|
| Gene vocabulary mismatch | High | Medium | Use intersection; 30 genes all protein-coding, likely in vocabulary |
| 239 samples too few for fine-tuning 316M params | Medium | High | Use V1-10M model; freeze early layers; heavy regularization |
| Human→bovine transfer fails | Medium | Medium | Focus on human results; note cross-species limitation |
| Geneformer installation issues | Low | Low | Docker image available; fallback to scGPT |

---

## 5. Recommended Priority

1. **Phase 1 first**: embedding extraction is fast and low-risk
2. **Phase 2 if Phase 1 succeeds**: full fine-tuning for the manuscript
3. **Phase 3 if Phase 2 succeeds**: cross-species is a strong contribution
4. **Phase 4 as stretch goal**: in silico perturbation is most novel

**Estimated total time:** 11-18 hours across all phases.

---

## 6. Key papers to cite

| Paper | How to Cite |
|-------|-------------|
| **Geneformer** (Theodoris et al., Nature 2023) | Foundation model methodology |
| **scGPT** (Cui et al., Nature Methods 2024) | Alternative foundation model comparison |
| **MEDL** (Andrade et al., 2024) | Batch effect deep learning validation |
| **TEDDY** (Chevalier et al., 2025) | Scaling laws for single-cell foundation models |
| **BMFM-RNA** (Dandala et al., 2025) | Benchmarking framework for transcriptomic FMs |
| **GeneMamba** (Qi et al., 2025) | Efficient alternative to transformer FMs |
