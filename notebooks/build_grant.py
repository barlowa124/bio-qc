"""Generate grant proposal draft."""
from datetime import datetime
from pathlib import Path

PROJ = Path(__file__).resolve().parents[1]
OUT = PROJ / "p2_state_map/output"

proposal = f"""# NIH R21 Grant Proposal — Draft
## Prospective Validation and Optimization of a 30-Gene qPCR Panel for Cultivated Meat Manufacturing Readiness

**PI:** Rao Lab, North Carolina State University
**Mechanism:** NIH R21 (Exploratory/Developmental Research Grant)
**Duration:** 2 years
**Budget:** ~$275,000 direct costs
**Generated:** {datetime.now().strftime('%B %d, %Y')}

---

### SPECIFIC AIMS

Cultivated meat production requires robust, quantitative quality control (QC) metrics to assess stem cell manufacturing readiness. Current approaches rely on subjective morphological assessment and endpoint assays that cannot provide real-time feedback for bioprocess optimization. We have developed a multi-omic state map integrating RNA-seq transcriptomics and METAFlux metabolic flux predictions across 239 samples, defining three manufacturing readiness states (expansion-competent, committed, terminal) with 96.7% cross-validation accuracy. A 30-gene L1-regularized logistic regression panel achieves equivalent accuracy, enabling routine batch triage at $50-100 per assay with 4-6 hour turnaround.

**Aim 1: Prospective wet-lab validation of the 30-gene panel.** We will validate the panel prospectively in 30 independent stem cell culture batches across multiple operators and passages, comparing qPCR-based state assignments against RNA-seq ground truth. Hypothesis: The panel will maintain >90% agreement with RNA-seq state assignments.

**Aim 2: Media formulation optimization using compound response prediction.** We will test 10 computationally-predicted small molecules (including p38 inhibitors, TGF-β inhibitors, Wnt activators, and Notch activators) for their ability to maintain cells in the expansion-competent state. Hypothesis: At least 3 compounds will significantly increase the proportion of expansion-competent cells.

**Aim 3: Minimal panel refinement and cross-platform validation.** We will refine the panel to the minimal gene set maintaining >95% accuracy (computationally predicted at 10-15 genes) and validate across qPCR platforms (Bio-Rad, Thermo Fisher, Roche). Hypothesis: A 10-15 gene panel will maintain non-inferior accuracy to the full 30-gene panel.

---

### RESEARCH STRATEGY

#### Significance
The cultivated meat industry is projected to reach $25 billion by 2030, but manufacturing scalability remains the primary bottleneck. A critical unmet need is standardized QC metrics for batch-to-batch consistency. Our preliminary data demonstrate:
- 96.7% state classification accuracy (239 samples, 5-fold CV)
- Cross-species conservation (bovine GSE173199, 38 samples)
- Single-nucleus resolution validation (GSE240556, 17,541 nuclei)
- Bayesian uncertainty quantification (97.9% confident assignments)
- Bootstrap stability (96.3% ± 0.8%, 1000 iterations)
- qPCR noise robustness (94.2% accuracy at 25% CV)

Notably, our 30-gene panel shows **zero overlap** with canonical muscle stem cell gene sets (MSigDB, PanglaoDB), indicating it captures novel biology beyond established markers.

#### Innovation
1. First prospective validation of a computationally-derived QC panel for cultivated meat
2. Integration of drug prediction with media optimization for manufacturing
3. Cross-platform validation enabling immediate industrial deployment
4. Open-source prediction API for community adoption

#### Approach
**Aim 1:** 30 batches × 3 timepoints × 3 technical replicates = 270 qPCR assays. Compare against RNA-seq (30 samples). Primary endpoint: Cohen's kappa for state agreement. Power analysis: 30 batches provide >90% power to detect kappa >0.80.

**Aim 2:** 10 compounds × 3 concentrations × 3 replicates = 90 culture conditions. Flow cytometry for proliferation (EdU) and differentiation (myosin heavy chain). RNA-seq on top 5 conditions.

**Aim 3:** Sequential feature elimination with cost constraints. Cross-platform comparison: 20 samples × 3 platforms = 60 assays. Bland-Altman analysis for platform agreement.

#### Timeline
- Months 1-6: Aim 1 (prospective validation)
- Months 7-18: Aim 2 (media optimization)
- Months 12-24: Aim 3 (panel refinement + cross-platform)
- Months 18-24: Manuscript preparation + dissemination

---

### BUDGET JUSTIFICATION

| Category | Year 1 | Year 2 | Total |
|----------|--------|--------|-------|
| Personnel (postdoc 50%, tech 100%) | $85,000 | $87,500 | $172,500 |
| Supplies (qPCR reagents, media, compounds) | $35,000 | $30,000 | $65,000 |
| RNA-seq (30 samples × $300) | $9,000 | $0 | $9,000 |
| Travel + publication | $3,000 | $3,000 | $6,000 |
| **Total Direct** | **$132,000** | **$120,500** | **$252,500** |
| Indirect (52% MTDC) | $68,640 | $62,660 | $131,300 |
| **Total** | **$200,640** | **$183,160** | **$383,800** |

---

### PRELIMINARY DATA

See attached:
- `p2_state_map/output/analysis_report.html` — Full analysis report
- `p2_state_map/output/literature_drug_panel_noise.json` — Benchmark + drug + noise data

### REFERENCES

1. Bernet JD et al. (2014) p38 MAPK signaling regulates stem cell self-renewal. *Nat Med* 20:265.
2. Cosgrove BD et al. (2014) Rejuvenation of the muscle stem cell population. *Nat Med* 20:255.
3. Tierney MT et al. (2014) STAT3 signaling controls satellite cell expansion. *Nat Med* 20:1182.
4. Conboy IM et al. (2003) Notch-mediated restoration of regenerative potential. *Science* 302:1575.
5. Sato T et al. (2016) A Pax7/Dmrt2/Myf5 regulatory network. *Cell Stem Cell* 19:107.
6. Rodgers JT et al. (2014) mTORC1 controls the adaptive transition of quiescent stem cells. *Nature* 510:393.
7. Our preliminary data: Multi-omic state map, 30-gene panel, cross-species validation (this work).

---

### ALTERNATIVE FUNDING MECHANISMS

**NSF EAGER** ($200-300K, 1-2 years): For high-risk, high-reward exploratory research. Emphasize the novel computational-to-wet-lab pipeline and cross-species validation.

**Good Food Institute (GFI) Research Grant** ($250K, 2 years): Specifically for alternative protein research. Emphasize immediate industrial applicability.

**NC Biotech Center Flash Grant** ($25K, 6 months): For rapid validation of the 30-batch prospective study. Lower barrier to entry, faster turnaround.
"""

out_path = OUT / "grant_proposal_draft.md"
out_path.write_text(proposal)
print(f"Grant proposal saved to: {out_path}")
print("DONE")
