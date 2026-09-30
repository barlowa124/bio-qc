# Installation

## Requirements

- Python 3.10+
- 16GB+ RAM recommended
- CUDA-capable GPU optional (for VAE training)

## Dependencies

```bash
pip install -r requirements.txt
```

Key packages:
- `numpy`, `pandas`, `scipy`: data processing
- `scikit-learn`: machine learning
- `torch`: deep learning (VAE)
- `anndata`, `scipy.sparse`: single-cell data
- `matplotlib`, `seaborn`: visualization

## Data

Public datasets are downloaded automatically or can be placed in:
- `cross_species_validation/`: bovine and porcine GEO data
- Source data paths configured in `shared_data/data_paths.yaml`
