# spatial_qc

Spot-level QC and filtering-strategy comparison for 10x Visium data.
Reads a Space Ranger export (`filtered_feature_bc_matrix/` MTX +
`spatial/tissue_positions_list.csv`), computes the standard per-spot
metrics (total counts, genes detected, mitochondrial fraction,
in-tissue), and benchmarks three spot-filtering strategies side by
side.

```bash
python -m spatialqc.cli SAMPLE_DIR --dataset NAME --out results/qc.json
```

## Strategies compared

| strategy | rule |
|---|---|
| `fixed` | static cutoffs: >=500 counts, >=200 genes, <=25% mito |
| `mad` | per-sample adaptive: drop spots >5 MADs below median on log counts/genes or above on mito (scverse-style outlier detection) |
| `tissue_only` | keep `in_tissue == 1`, ignore molecule evidence |

The comparison reports spots/counts retained, genes still detected,
median depth, and (when scanpy is installed) Leiden clusters and graph
modularity per strategy. Modularity is a structure-quality proxy, not
a ground-truth score.

## Committed result

`results/v1_adult_mouse_brain.json`: 10x Genomics V1 Adult Mouse Brain
public sample (2,702 spots, 32,285 genes, Space Ranger 1.1.0 export).

| strategy | spots kept | counts kept | clusters | modularity |
|---|---|---|---|---|
| fixed | 2,611 (96.6%) | 98.6% | 16 | 0.8368 |
| mad | 2,625 (97.2%) | 99.0% | 15 | 0.8406 |
| tissue_only | 2,702 (100%) | 100% | 14 | 0.8336 |

MAD filtering keeps more real spots than fixed cutoffs with
slightly better graph structure on this sample. The distribution-adaptive approach
pays on this sample. The `tissue_only` line is near-passthrough
because the *filtered* matrix already excludes out-of-tissue barcodes
upstream. Its discriminating power shows only on the unfiltered
(in-tissue + fiducial-frame) export, which is the same call for anyone
who has it.

## Limits

- One public sample. The strategy ranking is a demonstration of the
  measurement, not a claim about Visium data generally.
- Mito genes are detected by `mt-` symbol prefix, right for mouse and
  human nomenclature, wrong for organisms with different conventions.
- Modularity compares filtered subsets on their own graphs. It answers
  whether a filtering keeps structure-coherent spots, not which
  clusters are biologically right. There is no ground truth here.

## Tests

```bash
PYTHONPATH=src python -m pytest tests/ -q    # 7 tests, synthetic MTX fixture
```
