# statgen report

Cohort mode `demo`; trait `synthetic`. Synthetic and simulated-trait numbers exercise the pipeline and are not genetic findings.

## QC waterfall

800 -> 800 samples; 4000 -> 3840 variants.

| filter | axis | removed | kept |
|---|---|---:|---:|
| sample_missingness | sample | 0 | 800 |
| variant_missingness | variant | 0 | 4000 |
| min_maf | variant | 23 | 3977 |
| hwe | variant | 137 | 3840 |

## Stratification

10 PCs computed on 3840 variants; PC1 explains 6.0% of genotype variance and correlates with the ancestry label at |r| = 1.00.

## Association

linear test, 3840 variants on 800 samples with 4 covariates (2 PCs).

Genomic control $\lambda$ = 0.965 with covariates vs 2.646 without; 4 variants pass the Bonferroni threshold 1.30e-05.

Lead hit rs102113 (pos 2114): beta 0.443 (se 0.063), p = 3.54e-12.

Planted causal variants: 4 of 6 surviving QC recovered at Bonferroni (8 planted; 2 failed pooled HWE — the Wahlund effect of mixing groups; beta correlation 0.99.

Of the 4 Bonferroni hits, 0 are non-causal.

![manhattan](manhattan.png)

![qq](qq.png)

## Limits

Single-marker tests on a small cohort; no LD-aware fine-mapping or imputation. Simulated traits use recorded seeds and planted effects only to exercise recovery.
