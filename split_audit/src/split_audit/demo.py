"""Deterministic demo dataset with planted leakage.

The generator emits a fixed synthetic corpus: random scaffold/sequence
partitions plus a recorded set of planted cross-split duplicates. The
planted leaks are the ground truth the audit must recover, so the demo
is self-verifying — like statgen's planted-truth cohort, demo numbers
are mechanics checks, not findings.
"""
from __future__ import annotations

import random

_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"


def _rand_seq(rng: random.Random, length: int) -> str:
    return "".join(rng.choice(_ALPHABET) for _ in range(length))


def _mutate(rng: random.Random, seq: str, n: int) -> str:
    seq = list(seq)
    for pos in rng.sample(range(len(seq)), min(n, len(seq))):
        seq[pos] = rng.choice(_ALPHABET)
    return "".join(seq)


def make_demo(seed: int = 7, n_train: int = 40, n_test: int = 15,
              n_planted: int = 3) -> tuple[list[dict], list[dict], list[str]]:
    """Return (scaffold rows, sequence rows, planted sequence ids).

    Scaffolds: two scaffold ids are planted across the split boundary.
    Sequences: n_planted test sequences are near-duplicates (3 mutations)
    of train sequences; the rest are random and far below any k-mer
    threshold.
    """
    rng = random.Random(seed)

    scaf_rows = []
    for i in range(n_train + n_test):
        split = "train" if i < n_train else "test"
        scaf_rows.append({"id": f"mol-{i:03d}", "scaffold_id": f"scaf-{i}",
                          "split": split})
    planted_scafs = []
    for j, test_i in enumerate(rng.sample(range(n_train, n_train + n_test),
                                          2)):
        planted_scafs.append(scaf_rows[test_i]["scaffold_id"])
        scaf_rows[test_i]["scaffold_id"] = scaf_rows[j]["scaffold_id"]

    seq_rows = []
    train_seqs = [_rand_seq(rng, 60) for _ in range(n_train)]
    for i, seq in enumerate(train_seqs):
        seq_rows.append({"id": f"tr-{i:03d}", "sequence": seq,
                         "split": "train"})
    planted_ids = []
    donors = rng.sample(range(n_train), n_planted)
    for j in range(n_test):
        if j < n_planted:
            seq = _mutate(rng, train_seqs[donors[j]], 3)
            planted_ids.append(f"te-{j:03d}")
        else:
            seq = _rand_seq(rng, 60)
        seq_rows.append({"id": f"te-{j:03d}", "sequence": seq,
                         "split": "test"})
    return scaf_rows, seq_rows, planted_ids
