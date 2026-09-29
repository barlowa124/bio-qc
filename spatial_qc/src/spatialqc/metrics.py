"""Per-spot QC metrics over a loaded Visium sample.

Metrics follow the scverse/squidpy vocabulary: total_counts,
n_genes_by_counts, mito fraction (gene symbols beginning ``mt-``),
joined against in_tissue from the positions file. All arrays are
spot-aligned with the input matrix rows.
"""

from __future__ import annotations

import numpy as np


def spot_metrics(sample: dict) -> dict:
    m = sample["matrix"]
    counts = np.asarray(m.sum(axis=1)).ravel()
    n_genes = np.asarray((m > 0).sum(axis=1)).ravel()
    mito_cols = np.array(
        [g.lower().startswith("mt-") for g in sample["genes"]])
    if mito_cols.any():
        mito = np.asarray(m[:, mito_cols].sum(axis=1)).ravel()
    else:
        mito = np.zeros(m.shape[0])
    with np.errstate(invalid="ignore", divide="ignore"):
        mito_frac = np.where(counts > 0, mito / counts, 0.0)
    in_tissue = np.array(
        [int(sample["positions"].get(b, {}).get("in_tissue", 0))
         for b in sample["barcodes"]])
    return {"total_counts": counts, "n_genes": n_genes,
            "mito_frac": mito_frac, "in_tissue": in_tissue}


def summary(sample: dict, met: dict | None = None) -> dict:
    met = met or spot_metrics(sample)
    it = met["in_tissue"].astype(bool)
    out = {
        "n_spots": int(len(met["total_counts"])),
        "n_genes": int(len(sample["genes"])),
        "n_in_tissue": int(it.sum()),
        "frac_in_tissue": round(float(it.mean()), 4),
        "counts": {
            "total": int(met["total_counts"].sum()),
            "median_in_tissue": float(np.median(met["total_counts"][it]))
                if it.any() else 0.0,
            "median_out_tissue": float(np.median(met["total_counts"][~it]))
                if (~it).any() else 0.0,
        },
        "genes_per_spot_in_tissue": float(np.median(met["n_genes"][it]))
            if it.any() else 0.0,
        "median_mito_frac_in_tissue":
            float(np.median(met["mito_frac"][it])) if it.any() else 0.0,
    }
    return out
