# Cultivated meat multi-omic analysis pipeline

## A manufacturing-readiness state map with 30-gene qPCR QC panel

**barlowa124**. Developed independently while affiliated with the Rao Lab, North Carolina State University. The lab is not an author

---

### Key Results

| Metric | Value |
|--------|-------|
| Multi-omic state map accuracy | **96.7%** (5-fold CV) |
| 30-gene panel accuracy | **96.7%** ± 1.0% |
| Cross-species validation | Bovine 3-state conservation confirmed |
| snRNA-seq validation | 17,541 bovine nuclei, 21/30 panel genes detected |
| Bayesian confidence | 97.9% confident assignments |
| Assay cost | \$50-100/batch |
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
