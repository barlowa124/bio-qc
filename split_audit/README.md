# split_audit

Train/test split-leakage auditor for bio-ML datasets. The rest of the
portfolio asserts frozen or homology-separated splits. This package
checks them after the fact.

Two audit modes, both stdlib-only.

**scaffolds.** Exact scaffold-id audit. Rows supply a `scaffold_id`
(Bemis-Murcko via RDKit or any clustering). Any scaffold id present in
more than one split is reported along with the leaked items.

**sequences.** Cross-split pairs are scored by k-mer Jaccard
similarity and flagged above a configurable threshold. This is a cheap
identity proxy, not an alignment. Production-scale corpora want
mmseqs2 or cd-hit.

## Demo

`split-audit demo` generates a deterministic synthetic corpus with
planted leakage: two shared scaffold ids and three near-duplicate test
sequences. The committed run under `results/` recovers all three
planted sequences at 0.73 to 0.76 Jaccard against a 0.6 threshold, with no
false positives across 600 cross-split pairs. Demo numbers are a
mechanics check, not findings.

## Use

```sh
pip install -e .
split-audit demo
split-audit scaffolds data/mols.csv
split-audit sequences data/proteins.csv -k 3 --threshold 0.6
```

Inputs are CSV with `id`, `split`, and `scaffold_id`/`sequence` columns
(overridable via `--*-col` flags). Each run writes a JSON result and a
markdown report under `--out` (default `results/`).

## Scope

- Exact scaffold ids only. The auditor does not compute scaffolds
  itself, so it can verify a claimed split but not invent one.
- k-mer Jaccard approximates identity at sequence length. Short
  sequences and low-complexity regions can inflate it.
- Reports what crosses the boundary. Whether a flag means "leakage"
  versus "legitimately similar chemotypes" is a modeling decision the
  caller owns.
