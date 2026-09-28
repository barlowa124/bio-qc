# organoid-qc

Transcriptomic fidelity scoring for organoid cultures: how faithfully does a
scRNA-seq profile of an organoid reproduce the cell-type composition and
marker expression of the human tissue it is meant to model?

Standard organoid QC is manual and qualitative, a few marker stains and a
microscopy read. This pipeline makes it quantitative and reproducible:
cluster the organoid, correlate each cluster against reference tissue
cell-type centroids, and report **per-population** fidelity signals
(coverage for missing cell types, unmapped fraction for aberrant populations, and
marker detection for immature phenotypes) plus QC flags against config'd
thresholds.

## Quickstart

```bash
uv venv --python 3.11 && uv pip install -e ".[dev]"
.venv/bin/python -m snakemake --cores 2
```

Runs the deterministic **demo dataset** (synthetic organoid with known
ground truth: one missing cell type, one aberrant population) and writes:

- `results/fidelity.json` - per-cluster mappings, aggregate scores and flags
- `results/umap.png` - joint embedding, organoid cells colored mapped/unmapped

To score real data, set `dataset.mode: h5ad` in `config/config.yaml` and
point it at annotated `.h5ad` files (see `docs/data-sources.md` for
CELLxGENE/GEO access notes).

## Demo result

The demo culture intentionally fails QC. The pipeline
recovers the planted failure modes: `celltype_coverage` 0.667 vs the 0.75
threshold (stem-like population missing), `unmapped_fraction` 0.25 (one
aberrant cluster, correlation -0.02 to every reference centroid), stem
marker detection 0.0. `qc_pass: false`.

## Real-data result (`results/fidelity_intestine.json`)

`dataset.mode: fetch` ran the same pipeline on CELLxGENE Census data:
3,000 cells from the **human intestine organoid cell atlas** (353k-cell
dataset) vs 2,979 human small/large-intestinal tissue cells labeled
enterocyte / goblet cell / stem cell:

- `celltype_coverage` **1.0** - all three reference types hit by ≥1
  cluster. Most clusters map to stem cell (the atlas is
  stem/progenitor-dominated. Mapped-cluster correlations run
  0.36-0.86, median ~0.71, strong tail 0.74+)
- `unmapped_fraction` **0.075** - a few small clusters match nothing
  (best correlation ≤0.23)
- Marker detection: enterocyte 1.0, goblet 1.0, stem cell 0.667.
  One of three stem markers (LGR5/OLFM4/SMOC2) under-expressed relative
  to the tissue reference, consistent with immature organoid stem cells
- `qc_pass: true`

Interpretation stays conservative. This says the atlas sample's
transcriptome composition tracks intestinal epithelium, not that the
organoids function as intestine. See `docs/evaluation.md`.

## Interpretation

Fidelity ≠ function. A passing score means the transcriptome composition
matches the reference panel. It says nothing about barrier function,
secretion, or drug response. See `docs/evaluation.md` for the full list of
what a pass does not mean and known failure modes.

## Lab automation (Opentrons Flex)

`src/organoid_qc/automation/` contains a real Opentrons Flex protocol.
`flex_organoid_dosing.py` prepares an 8-point, 3-fold compound serial
dilution in a deep-well block and doses an organoid 96-well plate in
triplicate with vehicle-only control columns. `simulate.py` runs it
through the official Protocol Engine simulator, producing an audited
liquid-handling run log (278 steps: every aspirate/dispense/blow-out).
The same file runs on the physical robot. The tests assert step counts,
dose volumes, well mapping and that control columns stay untouched.
Requires `pip install -e .[automation]`.

## Structure

```
config/config.yaml      markers, thresholds, dataset mode, the single source of truth
workflow/Snakefile      demo_data/load -> preprocess -> score -> report
src/organoid_qc/
  data/synthetic.py     deterministic synthetic organoid + reference
  data/load.py          demo | h5ad modes
  qc/preprocess.py      filter -> normalize -> HVG -> PCA -> Leiden (shared gene space)
  qc/score.py           centroid correlation, cluster mapping, QC flags
  qc/report.py          joint UMAP + fidelity summary
  automation/           Opentrons Flex dosing protocol + simulation validation
tests/                  synthetic-data, scoring, and automation tests
docs/                   data sources, evaluation semantics
```

MIT licensed.

## Related work

- [cytof-qc](https://github.com/barlowa124/cytof-qc) applies the same AnnData/scanpy pipeline shape to mass-cytometry gates; its per-population agreement table is the same "report per cluster, not pooled" discipline this repo applies to fidelity scores.
