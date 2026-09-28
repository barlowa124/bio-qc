# Cross-Species Validation: Status & Next Steps

## Completed
- **Bovine count matrices downloaded**: GSE173199 (38 samples, 27,607 genes)
  - `sf_diff_counts.csv`: 18 samples, serum-free vs growth media
  - `timecourse_counts.csv`: 20 samples, D0-D7 differentiation timecourse
- **Metadata confirmed** for all 3 target datasets (bovine + porcine)
- **Analysis pipeline ready**: `notebooks/process_bovine.py`

## Blocked: Ensembl Gene Symbol Mapping
The bovine count matrix uses ENSBTAG IDs. Need to map to gene symbols for cross-species alignment. Ensembl REST API and FTP are returning errors from this network.

### Resolution options (try when network improves):
1. Run `notebooks/download_bovine_annot.py`: tries multiple Ensembl FTP URLs
2. Use `curl` to download: `curl -o bovine_genes.tsv.gz "https://ftp.ensembl.org/pub/release-99/tsv/bos_taurus/Bos_taurus.ARS-UCD1.2.99.gene.txt.gz"`
3. Use Python `requests` library instead of `urllib`
4. Use `biomart` package: `pip install biomart && python -m biomart`

### Once mapping is obtained:
```bash
python notebooks/process_bovine.py  # Will auto-detect the mapping file
```

## Prospective Validation
Full protocol at `cross_species_validation/validation_plan.md`. Ready to hand to wet lab team. No computational blockers.

## Files
```
cross_species_validation/
├── GSE173199/                  # Bovine count matrices (downloaded)
│   ├── sf_diff_counts.csv      # 18 samples
│   └── timecourse_counts.csv   # 20 samples
├── GSE240556_metadata.json     # Bovine snRNA-seq metadata
├── GSE206914_metadata.json     # Porcine scRNA-seq metadata
├── validation_plan.md          # Prospective lab protocol
└── STATUS.md                   # This file
```