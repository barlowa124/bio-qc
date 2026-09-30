# Cross-Species Validation & Prospective Lab Protocol

## Part A: Public Cross-Species Datasets

### Confirmed Available Datasets

#### Bovine (Bos taurus): directly relevant to cultivated meat

| Accession | Samples | Type | Context | Paper |
|-----------|---------|------|---------|-------|
| **GSE173199** | 38 | Bulk RNA-seq | Serum-free media for cultivated meat; satellite cell differentiation timecourse | Stout et al. 2021 |
| **GSE240556** | ~20 | snRNA-seq | Single-nucleus resolution of bovine satellite cell differentiation; identifies proliferating/differentiating/reserve cells | Nature Comms Biology 2024 |
| **GSE742077** | ~10 | scRNA-seq | Single-cell heterogeneity of cultured bovine satellite cells | Lyu et al. 2021, Frontiers in Genetics |

#### Porcine (Sus scrofa)

| Accession | Samples | Type | Context |
|-----------|---------|------|---------|
| **GSE206914** | ~15 | scRNA-seq + scATAC-seq | Integrative multi-omic myogenic differentiation | BMC Biology 2023 |
| **GSE162455** | ~12 | Bulk RNA-seq | Porcine muscle stem cell (pMuSC) differentiation day 0-4 | 2023 |
| *Nature Sci Data 2025* | 12 | Strand-specific RNA-seq | PSC proliferation (24h, 48h) + differentiation (18h, 28h) | 2025 |

### Download Protocol

```bash
# 1. Install prerequisites
pip install GEOparse pysradb

# 2. Download metadata
python -c "
import GEOparse
gse = GEOparse.get_GEO(geo='GSE173199', destdir='./data')
print(gse.metadata)
"

# 3. For count matrices, use SRA toolkit:
prefetch SRR14530716  # example SRR from GSE173199
fastq-dump --split-files SRR14530716

# 4. Quantify with Salmon (same pipeline as human data):
salmon quant -i bos_taurus_index -l A -1 sample_R1.fq -2 sample_R2.fq -o output
```

### Cross-Species Analysis Plan

1. **Download count matrices** for GSE173199 (bovine) and GSE206914 (porcine)
2. **Map to 1:1 orthologs** using Ensembl Compara (human↔bovine, human↔pig)
3. **Run METAFlux** on ortholog-mapped expression matrices
4. **Project onto human state map** using the existing PCA transformation
5. **Compare readiness state distributions** across species
6. **Validate 30-gene panel**: check if orthologs of the 30 human genes discriminate states in bovine/porcine

### Conserved Marker Genes for Cross-Species Validation

| Category | Genes |
|----------|-------|
| Pluripotency | POU5F1, SOX2, NANOG, MYC, KLF4 |
| Myogenesis | PAX7, MYF5, MYOD1, MYOG, MYH3, DES, MYL1, NEB |
| MSC markers | CD44, CD73 (NT5E), CD90 (THY1), CD105 (ENG) |
| Adipogenesis | PPARG, CEBPA, FABP4 |
| Housekeeping | GAPDH, ACTB, HPRT1 |

---

## Part B: Prospective Lab Validation Protocol

### Objective
Prospectively validate the 30-gene QC panel for batch triage in Rao Lab stem cell cultures.

### Study Design
- **Design:** Prospective cohort, blinded
- **Sample size:** 30 batches minimum
- **Duration:** 8-12 weeks
- **Operators:** 3 (inter-operator variability assessment)

### Sampling Schedule
Each batch sampled at 3 timepoints:
| Timepoint | Day | Purpose |
|-----------|-----|---------|
| T0 | Day 0 (seeding) | Baseline state |
| T1 | Day 3 (mid-expansion) | Proliferation checkpoint |
| T2 | Day 7 (pre-differentiation) | Final QC before differentiation commit |

### Measurements per Timepoint
1. **30-gene qPCR panel** (primary): ~$50-100/batch
2. **Morphology scoring** (current standard, for comparison)
3. **Viability** (trypan blue exclusion)
4. **Cell count** (for proliferation rate)
5. **Optional:** RNA-seq on a subset for full transcriptomic validation

### Endpoints
| Endpoint | Definition | Timing |
|----------|-----------|--------|
| **Primary:** Panel accuracy | Concordance between qPCR prediction and actual batch outcome | End of study |
| **Secondary:** Time-to-decision | Hours from sampling to go/no-go call | Per timepoint |
| **Secondary:** Cost per batch | Total QC cost vs batch value | End of study |
| **Secondary:** Operator agreement | ICC across 3 blinded operators | End of study |

### Ground Truth Definition
Batch outcome determined at day 14-21:
- **Pass:** >80% viability, >2 population doublings, successful differentiation
- **Fail:** <50% viability, <1 population doubling, or failed differentiation

### Analysis Plan
1. **ROC analysis**: qPCR panel prediction vs binary pass/fail outcome
2. **Confusion matrix**: 3-state prediction vs actual state
3. **Longitudinal trajectory**: state transitions T0→T1→T2 per batch
4. **Operator ICC**: intraclass correlation for panel scores across operators
5. **Cost-benefit**: (failed batch cost × early detection rate) vs (panel cost × n batches)

### Required Resources
| Item | Quantity | Est. Cost |
|------|----------|-----------|
| qPCR master mix | 30 batches × 3 timepoints × 30 genes | ~$1,500 |
| RNA extraction kits | 90 extractions | ~$900 |
| cDNA synthesis kits | 90 reactions | ~$600 |
| Primers (30 pairs) | 1 synthesis run | ~$500 |
| Culture consumables | 30 batches | ~$3,000 |
| **Total estimated:** | | **~$6,500** |

### qPCR Plate Layout (96-well)
```
Rows A-C: 30 target genes (triplicate)
Row D: Housekeeping genes (GAPDH, ACTB, HPRT1)
Row E: No-template controls
Row F: Inter-plate calibrator
```

### Expected Outcomes
- Panel accuracy >85% for pass/fail classification
- Time-to-decision <6 hours from sampling
- Operator ICC >0.8
- Cost savings: 1 early-detected failed batch saves ~$2,000-5,000 in wasted culture time

---

## Part C: Immediate Next Steps

1. **This week:** Download GSE173199 count matrix, run METAFlux, project onto state map
2. **This week:** Order 30 primer pairs for qPCR panel
3. **Week 1-2:** Pilot 5 batches through full protocol (validate workflow)
4. **Week 3-10:** Main study: 30 batches
5. **Week 11-12:** Analysis + manuscript revision with cross-species + prospective results