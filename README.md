# bio-qc

[![ci](https://github.com/barlowa124/bio-qc/actions/workflows/ci.yml/badge.svg)](https://github.com/barlowa124/bio-qc/actions/workflows/ci.yml)


Quality-control pipelines for single-cell, spatial, cytometry, and
multi-omic data. Each subdirectory is a self-contained package with its
own tests, config, and AGENTS.md.


## Where this sits in the portfolio

`bio-qc` is the **lab-data QC pipelines** repo: a dependency-free FCS parser (`fcs_io`), mass-cytometry QC in both Snakemake and Nextflow (`cytof_qc`, `nf`), spatial-transcriptomics QC (`spatial_qc`), and organoid fidelity scoring (`organoid_qc`). Sibling repos:
[trust-tools](https://github.com/barlowa124/trust-tools) (agent security
and evals), [bio-qc](https://github.com/barlowa124/bio-qc) (lab-data QC
pipelines), [lab-informatics](https://github.com/barlowa124/lab-informatics)
(lab data plumbing and integrity),
[llm-posttraining](https://github.com/barlowa124/llm-posttraining)
(training-stage behavior work),
[protein-ml](https://github.com/barlowa124/protein-ml) (protein fitness
ML), and [mol-ml](https://github.com/barlowa124/mol-ml) (small-molecule
ML).

## Packages

| Directory | What it does |
|---|---|
| `cytof_qc/` | Mass-cytometry QC and batch-alignment pipeline: arcsinh transform, drift checks, clustering, per-population metrics, a deployed review app, and a Snakemake DAG (`workflow/Snakefile`) that runs on generated FCS fixtures with no downloads. |
| `fcs_io/` | Dependency-free FCS 3.0/3.1 parser and writer with explicit vendor-quirk handling; wired into `cytof_qc` as the fallback FCS reader. |
| `nf/` | DSL2 Nextflow pipeline in nf-core module style: `FCSIO_DEMO` -> `FCSIO_PARSE` -> `FCS_STATS`, with `meta.yml`/`environment.yml`/`stub:` per module, nf-test coverage, and a committed `qc_summary.jsonl`. |
| `spatial_qc/` | Visium spot-level QC metrics plus a filtering-strategy benchmark (fixed cutoffs vs MAD-adaptive vs tissue-only), with a committed run on the public V1 Adult Mouse Brain export. |
| `organoid_qc/` | Organoid fidelity scoring: scRNA-seq organoid clusters vs tissue-reference centroids, per-cluster and per-cell-type fidelity, fail-closed QC flags. |
| `scrna_qc/` | scverse-based single-cell RNA QC: threshold filtering with an auditable waterfall, UMAP + Leiden, Wilcoxon markers, and a per-cluster QC table. Deterministic synthetic demo or public AnnData input (PBMC 3k / CELLxGENE-compatible), Snakemake DAG. |
| `statgen/` | Statistical genetics: genotype QC (missingness/MAF/HWE with a reconciling waterfall), stratification PCA, single-variant linear/logistic association with genomic-control lambda, and the same claims-check layer as scrna_qc binding every report number to a results JSON. Deterministic synthetic cohort (planted ancestry + causal variants) or 1000 Genomes chr22 real-data mode. |
| `cultivated_meat_multiomic/` | Multi-omic (RNA + metabolic flux) methods demo on public data: clustering, biomarker-panel selection with calibration and conformal intervals, ablation, drift monitoring, and cross-species checks on bovine/porcine muscle. A methods demonstration, not a manufacturing claim. |

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
