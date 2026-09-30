# statgen

Statistical-genetics DAG in the `scrna_qc` shape: config-driven
thresholds, a Snakemake pipeline, a deterministic synthetic cohort plus
a real public-data mode, committed result artifacts, and the same
claim-verification layer that audits every number in the markdown
report against recorded results.

```
data -> qc -> strat -> assoc -> report
```

## Run it

```bash
python -m pip install -e .[test]
snakemake -s workflow/Snakefile --cores 1          # demo cohort
PYTHONPATH=src python -m pytest tests/ -q
snakemake -n -s workflow/Snakefile                 # DAG dry run
```

## Demo cohort (synthetic, deterministic)

`dataset.mode: demo` generates 800 samples x 4,000 variants: two
ancestry groups differentiated by a Balding-Nichols model (Fst = 0.08),
2% missing calls, 8 planted causal variants at recorded effect sizes,
and a quantitative trait with an ancestry-shift confounder. The truth
lands in `results/planted_truth.json`, so recovery is auditable.

What the committed demo artifacts show:

- QC waterfall reconciles. 4,000 -> 3,840 variants (23 fail MAF, 137
  fail pooled HWE) and 1 of 800 samples fails the heterozygosity
  screen. Two of the planted causals die at HWE: mixing ancestry
  groups produces the Wahlund effect, and the filter removes exactly
  the stratified variants it should. Het outliers are z-scored within
  ancestry group, so a real population's het level is not mistaken for
  contamination.
- PC1 recovers the planted structure at |r| = 1.00 with the ancestry
  label on 6.0% of genotype variance.
- The relatedness scan (GRM off-diagonal on PC-projected genotypes,
  so ancestry sharing does not count as relatedness) flags 0 pairs
  above the 0.125 third-degree mark; a planted-duplicates test shows
  the same scan catching true replicates.
- PC covariates control inflation, dose-resolved: lambda_GC 2.63 with
  no PCs, 0.97 already at 1 PC, 0.97 at 2 PCs. One component absorbs
  the entire two-group confound.
- Signals come back without false positives. All 4 Bonferroni hits are
  planted causals (4 of 6 post-QC recovered. Beta correlation 0.99).
  The two misses are sub-threshold effects, reported as such.

Demo numbers exercise the pipeline. They are not genetic findings.

`config/demo_logistic.yaml` reruns the same cohort with the liability
thresholded to a binary trait and the Rao score test: committed under
`results/demo_logistic/`, lambda 1.72 -> 0.99 with PCs, 3 Bonferroni
hits, all causal (a thresholded trait trades power for the case/control
shape).

## Real-data mode (1000 Genomes chr22 head)

`config/kg_chr22.yaml` reads a 14 MB range-fetch of the phase3 chr22
call set (the first 20,000 biallelic variants across 2,504 real
samples) plus the integrated-call panel for superpopulation labels and
a real sex covariate. Fetch commands are in the config header. Nothing
public is committed beyond the derived artifacts under
`results/kg_chr22/`.

The call set ships no phenotype, so the config simulates a trait on the
real genotypes with a recorded seed. Causal variants are planted only
on common variants with limited continental frequency spread
(`max_anc_spread`), keeping the trait unconfounded. Setting
`ancestry_shift` > 0 deliberately reintroduces stratification.

Committed artifacts show what real data does to a GWAS pipeline:

- 20,000 -> 1,866 variants: 86% fail MAF >= 0.01 (the call set is dense
  with rare variation) and 971 fail pooled HWE across continental
  groups. 4 of 2,504 samples fail the within-superpopulation het screen.
- PC1 alone carries 4.5% of ancestry variance across five
  superpopulations. Real structure needs multiple PCs.
- lambda_GC by covariate PCs: 11.0 / 5.1 / 7.6 / 2.0 / 2.5 / 1.9 for
  0..5 PCs. Non-monotone because different PCs capture different
  structure than the planted trait loads on. It stays above 1 because a
  chr22-head slice is one LD-dense region and genome-wide lambda needs
  genome-wide sampling.
- The GRM relatedness scan flags 79,583 pairs at kinship > 0.125 with
  off-diagonal sd 0.058 (4x the demo's dispersion). 1,866 colocalized
  variants cannot estimate kinship reliably; the report prints the
  count with that caveat rather than letting it read as a cohort of
  cryptic relatives.
- 19 Bonferroni hits decompose as 2 recovered causals + 17 LD proxies
  within 250 kb + 0 unexplained. Two planted causals are perfectly
  collinear (identical statistics 619 bp apart). Their opposing true
  effects cancel marginally, which is the failure mode fine-mapping
  exists to solve.

## Claim verification

`src/statgen/claims.py` is the same verifier `scrna_qc` uses (ported
from `oncology_coscientist`), vendored byte-identical and sha256-pinned
by the repo-root parity test. This package is its second consumer.
Every number in `report.md` must re-derive from
`results/{qc_waterfall,strat_summary,assoc_summary,summary,planted_truth}.json`
at display tolerance. Unbound tokens fail the check and land in
`claims_check.json`.

## Layout

- `config/config.yaml`: demo cohort, thresholds, association settings
- `config/demo_logistic.yaml`: binary-trait score-test run
- `config/kg_chr22.yaml`: real-data mode with fetch commands
- `src/statgen/data.py`: cohort generation + VCF/panel loading
- `src/statgen/qc.py`: missingness/MAF/HWE filters, reconciling waterfall
- `src/statgen/strat.py`: genotype PCA, PC1-ancestry diagnostics
- `src/statgen/assoc.py`: linear (FWL) / logistic (score) scans, lambda_GC
- `src/statgen/report.py`: figures, summary, markdown, claim check
- `src/statgen/claims.py`: vendored claim verifier (second consumer)
- `workflow/Snakefile`: the DAG
- `results/` holds committed demo artifacts, `results/demo_logistic/`
  and `results/kg_chr22/` hold the binary-trait and real-data runs

## Limits

Single-marker marginal tests only, no LD-aware fine-mapping,
imputation, or mixed models. The logistic path is a score test, not a
fitted logistic regression. Simulated traits are pipeline exercises
with recorded seeds and planted effects, never findings.
