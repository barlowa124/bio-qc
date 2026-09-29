"""Visium loader: Space Ranger filtered MTX matrix + positions CSV.

Reads the standard export layout
(``filtered_feature_bc_matrix/{matrix.mtx,features.tsv,barcodes.tsv}.gz``
plus ``spatial/tissue_positions_list.csv``) into a spots x genes CSR
matrix aligned with per-spot positions. V1 positions files are
headerless; the six columns are barcode, in_tissue, array_row,
array_col, pxl_col, pxl_row.
"""

from __future__ import annotations

import csv
import gzip
from pathlib import Path

import numpy as np
import scipy.io
import scipy.sparse

POS_COLS = ["barcode", "in_tissue", "array_row", "array_col",
            "pxl_col", "pxl_row"]


def _lines(path: Path):
    op = gzip.open if path.suffix == ".gz" else open
    with op(path, "rt") as f:
        yield from f


def load_visium(sample_dir: str | Path) -> dict:
    """Return {'matrix': spots x genes CSR, 'genes': symbols,
    'gene_ids': ensembl, 'positions': dict per barcode}."""
    d = Path(sample_dir)
    mm_dir = d / "filtered_feature_bc_matrix"
    if not mm_dir.exists():
        mm_dir = d                       # allow pointing at the dir

    def find(name):
        for cand in (mm_dir / name, mm_dir / (name + ".gz")):
            if cand.exists():
                return cand
        raise FileNotFoundError(f"{name}[.gz] under {mm_dir}")

    barcodes = [l.strip() for l in _lines(find("barcodes.tsv"))]
    feats = [l.rstrip("\n").split("\t")
             for l in _lines(find("features.tsv"))]
    gene_ids = [f[0] for f in feats]
    genes = [f[1] if len(f) > 1 else f[0] for f in feats]

    mat_path = find("matrix.mtx")
    op = gzip.open if str(mat_path).endswith(".gz") else open
    with op(mat_path, "rb") as f:
        mat = scipy.io.mmread(f)         # genes x barcodes
    # MatrixMarket here may carry a header comment line; mmread handles it.
    matrix = scipy.sparse.csr_matrix(mat.T)   # -> spots x genes
    if matrix.shape[0] != len(barcodes):
        raise ValueError(
            f"barcode count {len(barcodes)} != matrix rows "
            f"{matrix.shape[0]}")

    positions = {}
    pos_path = d / "spatial" / "tissue_positions_list.csv"
    if pos_path.exists():
        for rec in csv.reader(_lines(pos_path)):
            if len(rec) >= 6:
                positions[rec[0]] = dict(zip(POS_COLS, rec))
    return {"matrix": matrix, "barcodes": barcodes, "genes": genes,
            "gene_ids": gene_ids, "positions": positions}
