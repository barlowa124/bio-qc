# bio-qc

[![ci](https://github.com/barlowa124/bio-qc/actions/workflows/ci.yml/badge.svg)](https://github.com/barlowa124/bio-qc/actions/workflows/ci.yml)


Quality-control pipelines for single-cell and cytometry data. Two related
projects merged into one repository, each a self-contained package with its
own tests, config, and commit history (imported via subtree merge).

## Packages

| Directory | What it does |
|---|---|
| `cytof_qc/` | Mass-cytometry QC and batch-alignment pipeline: arcsinh transform, drift checks, clustering, per-population metrics, a deployed review app, and a Snakemake DAG (`workflow/Snakefile`) that runs on generated FCS fixtures with no downloads. |
| `fcs_io/` | Dependency-free FCS 3.0/3.1 parser and writer with explicit vendor-quirk handling; wired into `cytof_qc` as the fallback FCS reader. |
| `spatial_qc/` | Visium spot-level QC metrics plus a filtering-strategy benchmark (fixed cutoffs vs MAD-adaptive vs tissue-only), with a committed run on the public V1 Adult Mouse Brain export. |
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
