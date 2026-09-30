# scrna_qc report — demo

- cells in: 900, genes in: 1200
- cells after filtering: 783, genes after filtering: 1199
- filter removals (cells / genes): min_genes_per_cell 72/0, max_pct_mt 45/0, min_cells_per_gene 0/1
- Leiden clusters: 4

## Per-cluster QC

| cluster | cells | % of cells | median counts | median genes | median %mt | top marker |
|---|---|---|---|---|---|---|
| 0 | 359 | 45.85% | 3024.0 | 840.0 | 1.81 | G0020 |
| 1 | 191 | 24.39% | 2910.0 | 839.0 | 1.84 | G0076 |
| 2 | 155 | 19.8% | 2953.0 | 831.0 | 1.82 | G0099 |
| 3 | 78 | 9.96% | 3102.5 | 846.0 | 1.82 | G0160 |

## Provenance

git 01360b8 · scanpy 1.11.5 · anndata 0.12.19
