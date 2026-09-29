# Evaluation: what the fidelity score does and does not show

## Signals reported (results/fidelity.json)

Per cluster: best-matching reference cell type, Pearson correlation to every
reference centroid, mapped/unmapped status.

Aggregate:

- **celltype_coverage**: fraction of reference types hit by ≥1 mapped
  cluster. Low coverage = missing populations (e.g. the culture never
  produces goblet cells).
- **unmapped_fraction**: fraction of organoid cells in clusters matching
  no reference type above `mapping_correlation_min`. Captures aberrant or
  off-target differentiation.
- **per_type_marker_detection**: fraction of each type's expected marker
  panel detected in its best-matching organoid cluster. Catches immature
  phenotypes that map by profile but lack markers.
- **qc_flags / qc_pass**: each signal vs its config threshold. The run
  passes only if all pass.

## What a pass does NOT mean

- Not functional equivalence. Transcriptome correlation says nothing about
  barrier function, secretion, morphology, or response to stimulus.
- Not batch-comparable across datasets. Scores depend on the reference
  composition, sequencing depth and preprocessing choices. Compare
  cultures scored against the *same* reference with the *same* config.
- Not a release criterion. This is a research QC signal for flagging
  cultures worth investigating, not a gate.

## Known failure modes

- Leiden resolution changes cluster granularity and therefore mapping.
  it is config-controlled, not tuned per-run.
- Correlation-to-centroid assumes reference types are transcriptionally
  distinct. Subtypes within one reference label blur into one centroid.
- Marker detection thresholds on mean expression are sensitive to dropout.
  a sparse cluster can under-report markers it expresses.
