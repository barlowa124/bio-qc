# scrna_qc report — pbmc3k

- cells in: 2700, genes in: 32738
- cells after filtering: 2643, genes after filtering: 13697
- filter removals (cells / genes): min_genes_per_cell 0/0, max_pct_mt 57/0, min_cells_per_gene 0/19041
- Leiden clusters: 7

## Per-cluster QC

| cluster | cells | % of cells | median counts | median genes | median %mt | top marker |
|---|---|---|---|---|---|---|
| 0 | 1172 | 44.34% | 2329.0 | 813.5 | 1.75 | RPS12 |
| 1 | 639 | 24.18% | 2617.0 | 944.0 | 2.29 | FTL |
| 2 | 433 | 16.38% | 1963.0 | 830.0 | 2.19 | NKG7 |
| 3 | 342 | 12.94% | 1764.0 | 675.5 | 2.06 | CD74 |
| 4 | 36 | 1.36% | 5297.0 | 1570.0 | 1.96 | HLA-DPA1 |
| 5 | 13 | 0.49% | 917.0 | 350.0 | 1.57 | PF4 |
| 6 | 8 | 0.3% | 10320.5 | 2702.0 | 1.82 | TYMS |

## Cluster-assignment confidence (Mondrian conformal)

- alpha: 0.1, calibration cells: 794, held-out cells: 1849
- pooled prediction-set coverage: 0.8831800973499189
- mean set size: 0.8853434288804759

| cluster | held-out cells | coverage | ci95 | qhat | set size | fallback |
|---|---|---|---|---|---|---|
| 0 | 820 | 0.8853658536585366 | [0.8617, 0.9054] | 0.006192505359649658 | 0.8902439024390244 | False |
| 1 | 447 | 0.8791946308724832 | [0.8457, 0.9062] | 0.001300215721130371 | 0.8791946308724832 | False |
| 2 | 303 | 0.900990099009901 | [0.8622, 0.9298] | 0.025535881519317627 | 0.900990099009901 | False |
| 3 | 239 | 0.9246861924686193 | [0.8841, 0.9518] | 0.005427539348602295 | 0.9246861924686193 | False |
| 4 | 25 | 0.48 | [0.3003, 0.665] | 0.005855262279510498 | 0.48 | True |
| 5 | 9 | 0.5555555555555556 | [0.2666, 0.8112] | 0.005855262279510498 | 0.5555555555555556 | True |
| 6 | 6 | 0.5 | [0.1876, 0.8124] | 0.005855262279510498 | 0.5 | True |

## Provenance

git 1dce091 · scanpy 1.11.5 · anndata 0.12.19
