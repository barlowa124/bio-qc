# scrna_qc report — pbmc3k

- cells in: 2700, genes in: 32738
- cells after filtering: 2643, genes after filtering: 13697
- filter removals (cells / genes): min_genes_per_cell 0/0, max_pct_mt 57/0, min_cells_per_gene 0/19041
- Leiden clusters: 8

## Per-cluster QC

| cluster | cells | % of cells | median counts | median genes | median %mt | top marker |
|---|---|---|---|---|---|---|
| 0 | 1175 | 44.46% | 2327.0 | 813.0 | 1.75 | RPS12 |
| 1 | 481 | 18.2% | 2331.0 | 860.0 | 2.23 | LYZ |
| 2 | 430 | 16.27% | 1964.0 | 831.0 | 2.19 | NKG7 |
| 3 | 342 | 12.94% | 1764.0 | 675.5 | 2.06 | CD74 |
| 4 | 158 | 5.98% | 3842.5 | 1268.5 | 2.41 | LST1 |
| 5 | 36 | 1.36% | 5297.0 | 1570.0 | 1.96 | HLA-DPA1 |
| 6 | 13 | 0.49% | 917.0 | 350.0 | 1.57 | PF4 |
| 7 | 8 | 0.3% | 10320.5 | 2702.0 | 1.82 | TYMS |

## Provenance

git 01360b8 · scanpy 1.11.5 · anndata 0.12.19
