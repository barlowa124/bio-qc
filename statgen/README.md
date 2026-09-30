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
  fail pooled HWE). Two of the planted causals die at HWE: mixing
  ancestry groups produces the Wahlund effect, and the filter removes
  exactly the stratified variants it should.
- PC1 recovers the planted structure at |r| = 1.00 with the ancestry
  label on 6.0% of genotype variance.
- PC covariates control inflation. Genomic-control lambda drops from
  2.65 (no covariates) to 0.97 with 2 PCs plus sex/age.
- Signals come back without false positives. All 4 Bonferroni hits are
  planted causals (4 of 6 post-QC recovered. Beta correlation 0.99).
  The two misses are sub-threshold effects, reported as such.

Demo numbers exercise the pipeline. They are not genetic findings.

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
  with rare variation) and 992 fail pooled HWE across continental
  groups.
- PC1 alone carries 4.5% of ancestry variance across five
  superpopulations. Real structure needs multiple PCs.
- lambda_GC 1.89 with 5 PC covariates, 11.66 without. A chr22-head
  slice is one LD-dense region, so residual inflation reflects
  correlated variants as much as stratification. Genome-wide lambda
  needs genome-wide sampling.
- 17 Bonferroni hits decompose as 2 recovered causals + 15 LD proxies
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
- `config/kg_chr22.yaml`: real-data mode with fetch commands
- `src/statgen/data.py`: cohort generation + VCF/panel loading
- `src/statgen/qc.py`: missingness/MAF/HWE filters, reconciling waterfall
- `src/statgen/strat.py`: genotype PCA, PC1-ancestry diagnostics
- `src/statgen/assoc.py`: linear (FWL) / logistic (score) scans, lambda_GC
- `src/statgen/report.py`: figures, summary, markdown, claim check
- `src/statgen/claims.py`: vendored claim verifier (second consumer)
- `workflow/Snakefile`: the DAG
- `results/` holds committed demo artifacts; `results/kg_chr22/` holds the real-data set

## Limits

Single-marker marginal tests only, no LD-aware fine-mapping,
imputation, or mixed models. The logistic path is a score test, not a
fitted logistic regression. Simulated traits are pipeline exercises
with recorded seeds and planted effects, never findings.
