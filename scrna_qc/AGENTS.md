# Project Guidance

- The demo generator is a deterministic synthetic dataset. It exists so
  the pipeline and tests run without a download. Never present demo
  numbers as biological findings.
- `mode: pbmc3k` loads the public PBMC 3k dataset via
  `sc.datasets.pbmc3k()` (cached by scanpy after first download). Never
  commit `.h5ad` files; commit only compact derived artifacts under
  `results/`.
- QC metrics are per-cell measurements, not biological release criteria.
  Report filtering effects per cluster, not only pooled. A pooled pass
  can hide a depleted population.
- `config/config.yaml` is the single source of truth for thresholds and
  embedding parameters; no hardcoded cutoffs in `src/`.
- Run `PYTHONPATH=src python -m pytest tests/ -q` and
  `snakemake -n -s workflow/Snakefile` after changes.
