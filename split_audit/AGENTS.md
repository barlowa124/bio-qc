# Project Guidance

- The demo generator is deterministic synthetic data with planted
  leakage. It exists so the audit and tests run without a download.
  Never present demo numbers as findings.
- The scaffold auditor checks scaffold ids only and does not compute
  them. The sequence auditor uses k-mer Jaccard, an identity proxy
  rather than an alignment. Do not describe either as an exhaustive
  leakage guarantee.
- `config/config.yaml` is the single source of truth for the k-mer
  length and Jaccard threshold, plus the demo sizes. No hardcoded
  cutoffs in `src/`.
- Run `PYTHONPATH=src python -m pytest tests/ -q` after changes. Commit
  refreshed `results/` artifacts whenever demo behavior changes.
