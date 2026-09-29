# bio-qc

Quality-control pipelines for single-cell and cytometry data. Two related
projects merged into one repository, each a self-contained package with its
own tests, config, and commit history (imported via subtree merge).

## Packages

| Directory | What it does |
|---|---|
| `cytof_qc/` | Mass-cytometry QC and batch-alignment pipeline: arcsinh transform, drift checks, clustering, per-population metrics, plus a deployed review app. |
| `organoid_qc/` | Organoid fidelity scoring: scRNA-seq organoid clusters vs tissue-reference centroids, per-cluster and per-cell-type fidelity, fail-closed QC flags. |

## Running tests

Each package is independent. From its directory:

```bash
cd cytof_qc && PYTHONPATH=src python -m pytest tests/ -q
```

Each subdirectory retains its own `AGENTS.md` with project-specific rules,
which still apply.

## Why one repo

Both answer the same question (does this population-level measurement
match its reference) over different measurement technologies, with the
same score-and-report shape and the same rule that pooled metrics must not
hide failed subpopulations.
