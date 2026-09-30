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

## Cluster-assignment confidence (Mondrian conformal)

- alpha: 0.1, calibration cells: 792, held-out cells: 1851
- pooled prediction-set coverage: 0.8427876823338736
- mean set size: 0.8465694219340897

| cluster | held-out cells | coverage | qhat | set size | fallback |
|---|---|---|---|---|---|
| 0 | 823 | 0.8578371810449574 | 0.0034123659133911133 | 0.8651275820170109 | False |
| 1 | 337 | 0.9347181008902077 | 0.01832479238510132 | 0.9347181008902077 | False |
| 2 | 301 | 0.8438538205980066 | 0.017805397510528564 | 0.8438538205980066 | False |
| 3 | 239 | 0.7364016736401674 | 0.002150416374206543 | 0.7364016736401674 | False |
| 4 | 111 | 0.7477477477477478 | 0.0527573823928833 | 0.7477477477477478 | False |
| 5 | 25 | 0.72 | 0.012701869010925293 | 0.76 | True |
| 6 | 9 | 0.4444444444444444 | 0.012701869010925293 | 0.4444444444444444 | True |
| 7 | 6 | 0.6666666666666666 | 0.012701869010925293 | 0.6666666666666666 | True |

## Provenance

git f600fa4 · scanpy 1.11.5 · anndata 0.12.19
