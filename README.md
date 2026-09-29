# cytof-qc

QC and unsupervised-clustering benchmark for mass cytometry (CyTOF) data,
run end to end on two public panels from the HDCytoData benchmark suite
(Levine et al. 2015, human bone marrow, manually gated populations):

- **Levine_13dim**: 13 markers, 24 gated populations, one donor
- **Levine_32dim**: 32 markers, 14 gated populations, two donors (H1/H2)
  acquired on different days, giving a real two-batch comparison

The pipeline parses raw FCS files, splits phenotypic markers from
QC-only channels (DNA, viability, acquisition bookkeeping), applies the
standard arcsinh transform, runs channel-, population-, and
acquisition-level QC, Leiden-clusters events through the scverse
toolchain (`anndata` + `scanpy`), and scores the clusters against the
published manual gates.

## Results

### Levine_13dim

- 167,044 events clustered (the full file, gated and ungated), 81,747 of
  them carrying manual-gate labels. Scoring runs only on the labeled
  subset. This matches the usual benchmark protocol: clustering gated
  cells only would inflate the agreement numbers.
- 24 gated populations, 21 Leiden clusters
- Agreement with manual gates: **ARI 0.867, NMI 0.867**
- Mature lineages map roughly one-to-one onto clusters. Gated progenitor
  populations partially merge, which is the expected failure mode.
  Progenitor gates overlap in marker space and manual gating uses
  hierarchy the clustering never sees.

### Levine_32dim (harder panel, two donors)

- 265,627 events clustered, 104,184 scored on gated events
- 14 gated populations, 39 Leiden clusters
- Agreement: **ARI 0.412, NMI 0.728**, lower than 13dim as expected on
  the progenitor-rich panel. Consistent across donors (H1 ARI 0.432,
  H2 ARI 0.491 after batch alignment).
- The measured per-channel batch shifts between donors are committed in
  the metrics file. The largest correction was DNA2 at 2.3 arcsinh
  units, a normalization/donor effect rather than a staining failure
- CD16+ NK cells, HSCs, and pDCs recover cleanly with recall above 0.98.
  Monocytes and CD4/CD8 T cells split across multiple clusters, the
  known hard case at this resolution

`results/<panel>/metrics.json` holds the full channel report, per-population
event counts, drift and batch reports, and the per-population
recall/purity table.

![UMAP: manual gates vs Leiden clusters](results/levine_13dim/umap_comparison.png)

## Run it

```bash
pip install -r requirements.txt
bash scripts/fetch_data.sh                            # both panels
python -m cytof_qc.pipeline --panel levine_13dim      # writes results/levine_13dim/
python -m cytof_qc.pipeline --panel levine_32dim      # writes results/levine_32dim/
python -m pytest tests/ -q                            # 35 tests
```

### Web service

Live: <https://cytof-qc-production.up.railway.app> — upload `.fcs` or
`.csv` events and get channel QC, drift flags, and Leiden clusters with
a UMAP view:

```bash
cd app && npm install && npm run build && cd ..    # build the React app once
uvicorn cytof_qc.service:app --port 8000           # http://127.0.0.1:8000
```

or `docker build -t cytof-qc . && docker run -p 8000:8000 cytof-qc`.
CI also publishes the image: `docker run -p 8000:8000
ghcr.io/barlowa124/cytof-qc:latest`.

Uploaded data has no manual gates, so the service reports descriptive
QC and cluster structure only — never agreement metrics it cannot
support. Requests are logged to `service_requests.jsonl` (path, status,
latency; no payload data).

## Layout

- `cytof_qc/io.py`: FCS parsing via `fcsparser`. Channel isotopes are
  renamed to marker names from `$PnN`/`$PnS` metadata. Two filename
  conventions are supported (`Marrow1_<pop>_cells.fcs` and
  `..._normalized_<pop>_<Hn>.fcs`). Residual files normalize to
  `NotGated`, and batch suffixes become a `batch` column.
- `cytof_qc/transform.py`: `arcsinh(x / cofactor)`, cofactor 5.
- `cytof_qc/qc.py`: channel-level negative/zero-event rates,
  per-population event-count flags, acquisition-order drift gate, and
  the phenotypic/QC-only channel split.
- `cytof_qc/compensation.py`: apply a measured spillover matrix to the
  event table (oxide/abundance crosstalk correction).
- `cytof_qc/batch.py`: per-channel between-batch median alignment with a
  report of every applied shift.
- `cytof_qc/cluster.py`: AnnData -> neighbors -> Leiden -> UMAP.
- `cytof_qc/mapping.py`: ARI/NMI, contingency table, Hungarian
  cluster-to-population matching, per-population recall/purity.
- `cytof_qc/pipeline.py`: the whole run and figure output.
- `cytof_qc/service.py` + `app/`: FastAPI upload endpoint and a
  React/TypeScript report viewer (canvas-rendered UMAP, channel QC and
  drift tables). Docker build bundles both.

## Honest scope

This is a research-grade benchmark on public reference data. It
demonstrates the standard analysis path (FCS -> transform -> QC ->
cluster -> compare to gates) but is not validated against a clinical or
production gating workflow, and the Leiden clustering is intentionally
unoptimized: no marker weighting and no per-population tuning.

Two scope limits on the QC gates: acquisition-drift flags on these
benchmark files measure file-to-file population composition (the files
are concatenated per population, not a continuous acquisition), so they
are reported with that caveat. Batch alignment is a per-channel
median shift, which centers batches but does not model nonlinear
warping.

## Related work

- [organoid-qc](https://github.com/barlowa124/organoid-qc) runs the same AnnData/scanpy path (transform, cluster, compare to reference labels) on organoid fidelity instead of gated CyTOF populations.
