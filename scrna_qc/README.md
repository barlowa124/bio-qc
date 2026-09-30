# scrna_qc

Single-cell RNA-seq QC pipeline built on scverse (`scanpy`/`anndata`).
Four stages (ingest, threshold filtering, embedding, reporting) wired as a
Snakemake DAG with a config-driven threshold set and a deterministic
synthetic demo so the whole thing runs without a download.

## What it produces

| Artifact | Contents |
|---|---|
| `results/filter_waterfall.json` | Per-step cell/gene counts in and out; retained + removed reconcile to inputs |
| `results/markers.csv` | Wilcoxon per-cluster marker scores |
| `results/cluster_qc.csv` | Per-cluster n_cells, fraction, median counts/genes/%mt, top marker |
| `results/qc_summary.json` | Summary + provenance (git sha, scanpy/anndata versions, config snapshot) |
| `results/report.md` | Human-readable report; every number binds to a recorded value |
| `results/claims_check.json` | Claim-check verdict for `report.md` (see below) |
| `results/cluster_conformal.json` | Mondrian conformal coverage of cluster-assignment confidence (see below) |
| `results/umap_set_size.png` | UMAP colored by prediction-set size on held-out cells (grey = calibration) |
| `results/umap.png` | UMAP colored by Leiden cluster |
| `data/processed/*.h5ad` | Filtered and embedded AnnData (gitignored) |

QC is reported per cluster, not only pooled. A pooled pass can hide a
depleted population.

## Verified claims

`report.md` goes through a claim check (`claims.py`, ported from the
oncology_coscientist draft verifier in trust-tools): every numeric token
in the markdown is re-derived against the flattened numeric leaves of
`filter_waterfall.json` + `qc_summary.json`, at the token's display
tolerance (a `%` token also binds to the fraction form). A number that
no recorded computation produced lands in `claims_check.json` under
`unbound_claims`. Cluster names and version strings are identifiers, not
claims, and are exempt. Every claim also records `bound_to` (the
artifact key it was derived from), so the verdict is an inspectable
binding map, not just a boolean. The committed verdicts: all numbers
bound on both runs (68 demo, 107 pbmc3k claims including the conformal
table).

## Cluster-assignment confidence

`cluster_confidence.py` turns the hard Leiden labels into calibrated
prediction sets. A logistic classifier on the PCA coordinates scores
each cell per cluster. The vendored `conformal.py` (the same shared copy
pinned in comp-tox, protstab, and cultivated-meat; this is its fourth
consumer) calibrates nonconformity per true cluster via `mondrian_qhat`,
with the pooled level as fallback for clusters under
`conformal.min_group` calibration cells.

Coverage is reported per cluster, not only pooled. The per-cluster
number is the one that matters, and each carries a Wilson 95% interval
so a 9-cell cluster reads differently from a 900-cell one. On the
committed runs the pooled coverage is 0.80 (demo) / 0.88 (pbmc3k) at
alpha 0.1, and the per-cluster table shows why pooled-only reporting
would mislead: large well-separated clusters sit near 0.86 to 0.93,
while the demo's smallest cluster falls to 0.22 on the global fallback
and pbmc3k's rare clusters run 0.44 to 0.72 on single-digit cell counts.
Pooled coverage below 1-alpha is expected under Mondrian mixing; the
guarantee lives at the per-group level for groups with their own qhat,
and single-split numbers carry realization noise on top. Low
conditional coverage means the classifier is rarely confident correctly
inside that cluster, which is itself a QC signal. `umap_set_size.png`
maps the prediction-set size onto the embedding, with held-out cells
colored and calibration cells grey.

One caveat the record now makes visible: the partition depends on which
Leiden backend ran, since leidenalg and igraph flavors partition
differently at the same seed and resolution. `qc_summary.json`
provenance records python, scanpy, anndata, leidenalg, and igraph
versions, so a shift like 8 to 7 clusters has a diagnosable cause
instead of silent drift.

## Committed example runs

- `results/`: the synthetic demo (900 cells, 4 planted populations,
  injected low-quality and high-mito cells).
- `results/pbmc3k/`: the real public 10x PBMC 3k dataset (2,700 cells,
  32,738 genes). The run filtered 57 high-mito cells and ~19k
  low-observation genes, then produced 7 Leiden clusters whose top
  markers land on canonical PBMC families (LYZ/S100A8 monocytes,
  NKG7/GZMA NK, CD74/HLA-DPA1 antigen-presenting, TYMS cycling). These
  clusters confirm the pipeline mechanics on real data. They are not a
  novel biological finding; the standard pbmc3k workflow produces the
  same grouping.

Run the pbmc3k config yourself (downloads ~5.5 MB once):

```bash
SCRNA_QC_CONFIG=config/pbmc3k.yaml snakemake -s workflow/Snakefile -c1
```

## Run

```bash
pip install -e .
snakemake -s workflow/Snakefile -c1   # or per-stage: python -m scrna_qc.data ...
```

## Config

`config/config.yaml` is the single source of truth.

- `paths`: `data_dir` (gitignored intermediates) and `results_dir`
  (committed artifacts) — `config/pbmc3k.yaml` overrides these so the
  real run lands under `results/pbmc3k/`.
- `dataset.mode`: `demo` (synthetic, deterministic, seeded), `pbmc3k`
  (public PBMC 3k set via `sc.datasets.pbmc3k()`, cached by scanpy), or
  `h5ad` (local path, compatible with CELLxGENE/GEO/HCA exports).
- `qc`: `min_genes_per_cell`, `max_pct_mt`, `min_cells_per_gene`,
  `mt_prefix`.
- `embed`: normalization target, HVG count, PCs, neighbors, Leiden
  resolution, markers per cluster, seed.

Override the config path with `SCRNA_QC_CONFIG`.

## Tests

```bash
PYTHONPATH=src python -m pytest tests/ -q   # 35 tests
snakemake -n -s workflow/Snakefile          # DAG dry-run
```

## Scope

- QC metrics and clusters here are pipeline outputs, not biological
  findings. The demo dataset is synthetic; its "types" are sampling
  artifacts with marker blocks, not cell identities.
- Filtering removes cells and genes. `filter_waterfall.json` exists so
  the removal stays visible and auditable.
