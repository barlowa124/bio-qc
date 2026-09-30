# statgen report

Cohort mode `vcf`; trait `real_genotypes_simulated_trait`. Synthetic and simulated-trait numbers exercise the pipeline and are not genetic findings.

## QC waterfall

2504 -> 2500 samples; 20000 -> 1899 variants.

| filter | axis | removed | kept |
|---|---|---:|---:|
| sample_missingness | sample | 0 | 2504 |
| heterozygosity | sample | 4 | 2500 |
| variant_missingness | variant | 0 | 20000 |
| min_maf | variant | 17130 | 2870 |
| hwe | variant | 971 | 1899 |

## Stratification

10 PCs computed on 1899 variants; PC1 explains 9.0% of genotype variance

Relatedness scan (GRM off-diagonal after projecting out the 10 PCs, so ancestry sharing does not count as relatedness): max 1.150, sd 0.058, 79583 pairs above the 0.125 flag — with few regional variants the estimate is overdispersed; a high count reflects LD, not cryptic relatedness.; 5 ancestry groups capture 4.6% of its variance (eta^2).

## Association

linear test, 1899 variants on 2500 samples with 6 covariates (5 PCs).

Genomic control $\lambda$ = 1.856 with covariates vs 11.019 without; 19 variants pass the Bonferroni threshold 2.63e-05.

$\lambda$ by covariate PCs: 0 PCs -> 11.037, 1 PCs -> 5.050, 2 PCs -> 7.583, 3 PCs -> 2.007, 4 PCs -> 2.451, 5 PCs -> 1.856. Effect sizes are per copy of the counted allele (a1).

Lead hit 22:17003990 (pos 17003990): beta -0.466 (se 0.057), p = 7.33e-16.

Planted causal variants: 2 of 5 surviving QC recovered at Bonferroni (6 planted; 1 failed pooled HWE — the Wahlund effect of mixing groups; beta correlation 0.64.

Of the 19 Bonferroni hits, 17 sit within 250 kb of a planted causal (LD proxies) and 0 are unexplained — marginal tests cannot separate correlated variants, so per-variant attribution needs fine-mapping.

![manhattan](manhattan.png)

![qq](qq.png)

## Limits

Single-marker tests on a small cohort; no LD-aware fine-mapping or imputation. Simulated traits use recorded seeds and planted effects only to exercise recovery.
