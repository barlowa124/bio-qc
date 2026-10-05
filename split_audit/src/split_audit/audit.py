"""Split-leakage auditing for bio-ML datasets.

Two audit modes, both stdlib-only:

- ``audit_scaffolds``: molecules sharing a scaffold id must not span a
  train/test boundary. Exact-id audit; the caller supplies scaffold ids
  (Bemis-Murcko via RDKit, or any clustering) so the auditor stays
  dependency-free.
- ``audit_sequences``: cross-split pairs are scored by k-mer Jaccard
  similarity and flagged above a threshold. k-mer Jaccard is a cheap
  identity proxy, not an alignment; pair it with mmseqs2/cd-hit for
  production-scale corpora.

Both modes return a plain dict suitable for JSON serialization so the
report layer can bind every printed number to this file's output.
"""
from __future__ import annotations

import itertools
from typing import Iterable


def audit_scaffolds(rows: Iterable[dict], id_col: str = "id",
                    scaffold_col: str = "scaffold_id",
                    split_col: str = "split") -> dict:
    """Report scaffold ids that appear in more than one split."""
    by_scaffold: dict[str, dict] = {}
    n_rows = 0
    splits: set[str] = set()
    for row in rows:
        n_rows += 1
        split = row[split_col]
        splits.add(split)
        bucket = by_scaffold.setdefault(row[scaffold_col],
                                        {"splits": set(), "ids": []})
        bucket["splits"].add(split)
        bucket["ids"].append(row[id_col])

    shared = {scaf: b for scaf, b in by_scaffold.items() if len(b["splits"]) > 1}
    leaked_items = sorted({i for b in shared.values() for i in b["ids"]})
    return {
        "n_rows": n_rows,
        "n_scaffolds": len(by_scaffold),
        "n_shared_scaffolds": len(shared),
        "n_leaked_items": len(leaked_items),
        "leaked_items": leaked_items,
        "shared_scaffolds": sorted(shared),
        "scaffold_disjoint": not shared,
    }


def _kmers(seq: str, k: int) -> set[str]:
    return {seq[i:i + k] for i in range(len(seq) - k + 1)}


def audit_sequences(rows: Iterable[dict], id_col: str = "id",
                    seq_col: str = "sequence", split_col: str = "split",
                    k: int = 3, threshold: float = 0.6,
                    report_top: int = 10) -> dict:
    """Flag cross-split sequence pairs at or above a k-mer Jaccard cutoff."""
    by_split: dict[str, list[tuple[str, set[str]]]] = {}
    n_rows = 0
    for row in rows:
        n_rows += 1
        kmers = _kmers(row[seq_col].upper(), k)
        by_split.setdefault(row[split_col], []).append((row[id_col], kmers))

    splits = sorted(by_split)
    flagged, max_per_test = [], {}
    n_pairs = 0
    for a, b in itertools.combinations(splits, 2):
        for id_a, kmers_a in by_split[a]:
            for id_b, kmers_b in by_split[b]:
                n_pairs += 1
                union = len(kmers_a | kmers_b)
                sim = len(kmers_a & kmers_b) / union if union else 0.0
                for split, seq_id in ((a, id_a), (b, id_b)):
                    max_per_test[seq_id] = max(sim, max_per_test.get(seq_id, 0.0))
                if sim >= threshold:
                    flagged.append({"a": id_a, "b": id_b, "similarity": sim})

    flagged.sort(key=lambda p: -p["similarity"])
    flagged_items = sorted({p["a"] for p in flagged}
                           | {p["b"] for p in flagged})
    return {
        "n_rows": n_rows,
        "splits": splits,
        "k": k,
        "jaccard_threshold": threshold,
        "n_cross_split_pairs": n_pairs,
        "n_flagged_pairs": len(flagged),
        "flagged_pairs": flagged[:report_top],
        "n_flagged_items": len(flagged_items),
        "flagged_items": flagged_items,
        "leakage_free": not flagged,
    }
