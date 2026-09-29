# nf - Nextflow FCS QC pipeline

A DSL2 pipeline in nf-core module style that parses flow-cytometry FCS
files and emits per-channel QC summaries.

```
FCSIO_DEMO -> FCSIO_PARSE -> FCS_STATS
```

- `FCSIO_DEMO` writes deterministic gated-population FCS fixtures
  (Levine-style `<sample>_<population>_cells.fcs`) via the `fcsio`
  writer, so the DAG runs with no download.
- `FCSIO_PARSE` runs one task per file: FCS 3.x to `events.tsv` and
  `meta.json`, using the sibling `fcs_io` package (validated bit-exact
  against `fcsparser` on the Levine benchmark. See
  `fcs_io/validation/`).
- `FCS_STATS` computes per-channel mean/median/p95/negative-fraction
  JSON per sample. The entrypoint collects them into
  `results/qc_summary.jsonl` and then writes
  `results/run_manifest.json`: git sha, artifact sha256s, params,
  and tool versions in the `bio-qc/run-manifest@1` schema shared with
  the Snakemake DAG (`../scripts/run_manifest.py`).

## Run

```bash
nextflow run nf/main.nf                       # demo fixtures
nextflow run nf/main.nf --input /path/to/fcs  # real files, same DAG
nextflow run nf/main.nf -stub-run             # DAG structure only
```

`--demo-events N` sets events per population (default 200).

## Tests

```bash
nf-test test nf/tests/main.nf.test \
    nf/modules/local/fcsio_parse/tests/main.nf.test
```

`tests/main.nf.test` runs the full workflow (13 tasks) and checks the
emitted stats. The module test exercises `FCSIO_PARSE` on a committed
fixture plus a `-stub` run.

## Conventions

Modules carry `meta.yml`, `environment.yml`, `versions.yml` output, and
a `stub:` section per the nf-core module spec. `fcsio` is not on PyPI or
conda, so processes export `PYTHONPATH` pointing at `../fcs_io/src`
rather than declaring a package dependency. The `container` directive
is a plain `python:3.12-slim` for environments that run containers.

## Limitations

- Demo fixtures are deterministic and perfectly gated. They validate
  the DAG, not biology.
- Vendor filenames with spaces are handled via `stageAs`. Sample ids
  are sanitized to underscores in the summary.
- The summary covers channel distributions only. Cluster-level QC
  (arcsinh transform, clustering, ARI) stays in the Snakemake path under
  `cytof_qc/workflow/`, which carries the heavier scverse dependencies.
