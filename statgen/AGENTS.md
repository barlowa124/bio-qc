# Project Guidance

- The demo cohort is deterministic synthetic data (two ancestry groups,
  planted causal variants at recorded effect sizes in
  `results/planted_truth.json`). It exists so the pipeline and tests
  run without a download. Never present demo numbers as genetic
  findings.
- `mode: vcf` reads real genotypes (1000 Genomes chr22 head region by
  default) and simulates a trait on them with a recorded seed. The call
  set ships no phenotype. The simulated trait is stated in the report.
  Fetch commands live in `config/kg_chr22.yaml`. Never commit the VCF
  or panel files.
- `config/config.yaml` is the single source of truth for thresholds
  (missingness, MAF, HWE) and association settings. No hardcoded
  cutoffs in `src/`.
- `src/statgen/claims.py` is vendored byte-identical with
  `scrna_qc/src/scrna_qc/claims.py`. The root parity test pins the
  sha256. Edit both copies together.
- Association stats are single-marker marginal tests. LD between
  variants inflates hit counts and blurs causal attribution. Report
  lambda_GC with and without covariates, and don't claim more than the
  numbers show.
- Run `PYTHONPATH=src python -m pytest tests/ -q` and
  `snakemake -n -s workflow/Snakefile` after changes. Snakemake does
  not track source edits, so use `-F` to force reruns.
