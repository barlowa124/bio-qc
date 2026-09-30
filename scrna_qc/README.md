# scrna_qc

Single-cell RNA-seq QC pipeline built on scverse (`scanpy`/`anndata`).
Four stages (ingest, threshold filtering, embedding, reporting) run as a
Snakemake DAG with a config-driven threshold set and a deterministic
synthetic demo so the whole thing runs without a download.

## What it produces

| Artifact | Contents |
|---|---|
| `results/filter_waterfall.json` | Per-step cell/gene counts in and out; retained + removed reconcile to inputs |
| `results/markers.csv` | Wilcoxon per-cluster marker scores |
| `results/cluster_qc.csv` | Per-cluster n_cells, fraction, median counts/genes/%mt, top marker |
| `results/qc_summary.json` | Summary + provenance (git sha, scanpy/anndata versions, config snapshot) |
| `results/umap.png` | UMAP colored by Leiden cluster |
| `data/processed/*.h5ad` | Filtered and embedded AnnData (gitignored) |

QC is reported per cluster, not only pooled. A pooled pass can hide a
depleted population.

## Run

```bash
pip install -e .
snakemake -s workflow/Snakefile -c1   # or per-stage: python -m scrna_qc.data ...
```

## Config

`config/config.yaml` is the single source of truth.

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
PYTHONPATH=src python -m pytest tests/ -q   # 17 tests
snakemake -n -s workflow/Snakefile          # DAG dry-run
```

## Scope

- QC metrics and clusters here are pipeline outputs, not biological
  findings. The demo dataset is synthetic; its "types" are sampling
  artifacts with marker blocks, not cell identities.
- Filtering removes cells and genes. `filter_waterfall.json` exists so
  the removal stays visible and auditable.
