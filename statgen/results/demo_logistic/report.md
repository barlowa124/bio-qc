# statgen report

Cohort mode `demo`; trait `synthetic`. Synthetic and simulated-trait numbers exercise the pipeline and are not genetic findings.

## QC waterfall

800 -> 799 samples; 4000 -> 3840 variants.

| filter | axis | removed | kept |
|---|---|---:|---:|
| sample_missingness | sample | 0 | 800 |
| heterozygosity | sample | 1 | 799 |
| variant_missingness | variant | 0 | 4000 |
| min_maf | variant | 23 | 3977 |
| hwe | variant | 137 | 3840 |

## Stratification

10 PCs computed on 3840 variants; PC1 explains 6.0% of genotype variance

Relatedness scan (GRM off-diagonal after projecting out the 10 PCs, so ancestry sharing does not count as relatedness): max 0.079, sd 0.015, 0 pairs above the 0.125 flag. and correlates with the ancestry label at |r| = 1.00.

## Association

logistic test, 3840 variants on 799 samples with 4 covariates (2 PCs).

Genomic control $\lambda$ = 0.986 with covariates vs 1.721 without; 3 variants pass the Bonferroni threshold 1.30e-05.

$\lambda$ by covariate PCs: 0 PCs -> 1.731, 1 PCs -> 0.988, 2 PCs -> 0.986. Effect sizes are per copy of the counted allele (a1).

Lead hit rs101142 (pos 1143): beta -0.634 (se 0.110), p = 8.50e-09.

Planted causal variants: 3 of 6 surviving QC recovered at Bonferroni (8 planted; 2 failed pooled HWE — the Wahlund effect of mixing groups; beta correlation 0.97.

Of the 3 Bonferroni hits, 0 are non-causal.

![manhattan](manhattan.png)

![qq](qq.png)

## Limits

Single-marker tests on a small cohort; no LD-aware fine-mapping or imputation. Simulated traits use recorded seeds and planted effects only to exercise recovery.
