"""Spot-filtering strategies — the axis the benchmark compares.

Each strategy returns a boolean keep-mask over spots:

- ``fixed``: classic static cutoffs (min counts, min genes, max mito
  fraction). Universally used, blind to the sample's own distribution.
- ``mad``: per-sample adaptive thresholds — spots more than ``nmad``
  MADs below the median on log(counts)/log(genes), or above on mito
  fraction, are dropped. This is the scverse-recommended outlier
  detection and adapts to platform depth.
- ``tissue_only``: keep ``in_tissue == 1`` and nothing else — the
  fiducial-frame answer, ignoring molecule evidence entirely.
"""

from __future__ import annotations

import numpy as np

# Visium-scale defaults; not universal — see README for scope.
FIXED_MIN_COUNTS = 500
FIXED_MIN_GENES = 200
FIXED_MAX_MITO = 0.25
MAD_NMAD = 5


def fixed(met: dict, min_counts: int = FIXED_MIN_COUNTS,
          min_genes: int = FIXED_MIN_GENES,
          max_mito: float = FIXED_MAX_MITO) -> np.ndarray:
    return ((met["total_counts"] >= min_counts)
            & (met["n_genes"] >= min_genes)
            & (met["mito_frac"] <= max_mito))


def _mad_outlier(x: np.ndarray, nmad: int, side: str) -> np.ndarray:
    med = np.median(x)
    mad = np.median(np.abs(x - med)) or 1e-9
    d = np.abs(x - med) / mad
    if side == "lower":
        return (x < med) & (d > nmad)
    if side == "upper":
        return (x > med) & (d > nmad)
    return d > nmad


def mad(met: dict, nmad: int = MAD_NMAD) -> np.ndarray:
    """Drop spots that are distribution outliers on any of the three
    metrics; log1p on count axes per scverse convention."""
    lc = np.log1p(met["total_counts"])
    lg = np.log1p(met["n_genes"])
    bad = (_mad_outlier(lc, nmad, "lower")
           | _mad_outlier(lg, nmad, "lower")
           | _mad_outlier(met["mito_frac"], nmad, "upper"))
    return ~bad


def tissue_only(met: dict) -> np.ndarray:
    return met["in_tissue"].astype(bool)


STRATEGIES = {"fixed": fixed, "mad": mad, "tissue_only": tissue_only}
