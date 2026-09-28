# Data sources

## Demo mode (default)

`dataset.mode: demo` generates a deterministic synthetic pair: a reference
"tissue" with one population per marker key in `config.markers`, and an
organoid culture that recovers all-but-one type at reduced marker expression
(immature phenotype) plus one aberrant population matching no reference type.

The demo exists so the pipeline and tests run without downloads. Its scores
are fixture values, never biological findings.

## Real data: fetch mode (CELLxGENE Census)

`dataset.mode: fetch` pulls directly from CELLxGENE Census via
`organoid_qc.data.fetch`: a subsample of the human intestine organoid cell
atlas (dataset `776a1e4a…`, 353k cells) and a reference slice of human
small/large-intestinal tissue filtered to the panel's cell types. Sample
sizes, census version and filters are under `fetch:` in config. The run
is reproducible modulo census version (pinned).

Notes from the real run:

- Census obs filters scan the full human obs table (~75M cells). The data
  stage takes ~25 min. Subsequent stages take seconds.
- Cell-type labels are census-normalized, not the atlas's own labels.
  For example, intestinal stem cells are `stem cell`, not `intestinal stem cell`.
  Check `reference.obs.cell_type.value_counts()` before writing a panel.
- Census `var` indexes by Ensembl `feature_id`. `fetch.py` renames vars to
  HGNC `feature_name` so config marker panels resolve.

## Real data: h5ad mode

Set `dataset.mode: h5ad` and provide two annotated AnnData files:

- `reference_h5ad`: human tissue scRNA-seq with a cell-type column
  (default `cell_type`; set `obs_columns.cell_type`).
- `organoid_h5ad`: organoid scRNA-seq, Leiden clusters are computed in
  preprocessing, or supply your own via `obs_columns.cluster`.

Practical sources:

- **CELLxGENE Census** (`cellxgene-census` package): query human tissue by
  `tissue`/`cell_type`, export `.h5ad` slices. Curated, consistently
  annotated, the best reference source.
- **GEO**: organoid atlas papers deposit count matrices + metadata. Assemble
  into AnnData. Check each study's own cell-type labels before merging.
- **Human Cell Atlas data portal**: larger reference collections, more
  assembly work.

Update `config.markers` to the tissue's canonical marker panel. The panel
is part of the experiment definition, not a default to reuse blindly.
