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

## Run receipts

Every run also emits `results/run_receipts.jsonl` — a hash-chained
execution record in the same format as `trust-tools/inference_receipts`:
one record per task binding the task's staged inputs (symlink targets,
sha256), written outputs (sha256), an env fingerprint (the task's
`versions.yml` plus the engine string), the Nextflow task hash, and
start/complete timestamps. `chain_prev` links records so a deleted,
reordered, or tampered task breaks the chain; verify with:

```bash
python ../scripts/run_receipts.py --verify results/run_receipts.jsonl
```

Receipts are written by `workflow.onComplete` (the trace writer flushes
before the handler runs) and can be regenerated deterministically after
any run:

```bash
python ../scripts/run_receipts.py --trace results/trace.txt \
    --out results/run_receipts.jsonl --run-name <name>
```

The receipts are a durable record of what ran. Two verification levels
exist: `--verify` checks record integrity (chain + embedded hashes),
while `--rehash` re-computes recorded files where they still exist —
inputs resolve via the staged symlink targets, outputs via
`workdir`/name — and reports hash mismatches or files cleaned away.
Receipts are emitted by the entry workflow's `onComplete`, so the named
`FCS_QC` workflow exercised by nf-test does not produce them; the
committed artifact comes from a real `nextflow run main.nf`.

## Run

```bash
nextflow run nf/main.nf                       # demo fixtures
nextflow run nf/main.nf --input /path/to/fcs  # real files, same DAG
nextflow run nf/main.nf -stub-run             # DAG structure only
```

`--demo-events N` sets events per population (default 200).

## Tests

```bash
cd nf   # module tests resolve bin/ from the project dir — run from nf/
nf-test test tests/main.nf.test \
    modules/local/fcsio_parse/tests/main.nf.test \
    modules/local/fcs_stats/tests/main.nf.test \
    modules/local/fcsio_demo/tests/main.nf.test
```

`tests/main.nf.test` runs the full workflow (13 tasks) and checks the
emitted stats. Each module has its own process-level test plus a
`-stub` run. `versions.yml` assertions check structure rather than
md5, since the recorded python version varies by host.

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
