# Project Guidance

- Use only public datasets (CELLxGENE Census, GEO, Human Cell Atlas) or the
  built-in deterministic synthetic demo. Never commit raw downloaded archives
  or `.h5ad` files; commit only compact derived artifacts under `results/`.
- Do not claim biological validity of any fidelity score beyond what the data
  supports. Scores are research-grade QC signals, not release criteria and
  not evidence of organoid function.
- The synthetic demo exists so the pipeline and tests run without large
  downloads; never present demo numbers as biological findings.
- Fidelity scores must be reported per-cluster and per-cell-type, not only as
  one pooled number — pooled scores hide missing and aberrant populations.
- `config/config.yaml` is the single source of truth for reference markers,
  scoring weights, and QC thresholds; no hardcoded biology in `src/`.
- Run `python -m pytest tests/` and `snakemake -n` after changes.
