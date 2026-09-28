"""Spillover compensation for mass-cytometry event tables.

Mass cytometry has far less channel overlap than fluorescence, but isotope
crosstalk still exists: oxide formation leaks signal into the M+16 channel,
and abundance sensitivity leaks into M-1/M+1. A spillover matrix ``S``
describes the fraction of each source channel's signal observed in each
target channel; compensation solves ``compensated = observed @ inv(S)``.

The matrix must be measured from single-stain controls by the caller; this
module only applies it and refuses matrices that are not invertible or do
not match the channel list.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def compensate(
    events: pd.DataFrame,
    channels: list[str],
    spillover: np.ndarray,
) -> pd.DataFrame:
    """Apply compensation and return a copy of `events` corrected.

    ``spillover[i, j]`` is the fraction of channel ``i``'s true signal
    observed in channel ``j`` (diagonal entries are 1.0). Rows/columns are
    ordered as ``channels``. Non-square or misordered input is rejected;
    the inverse must exist and produce finite values.
    """
    S = np.asarray(spillover, dtype=float)
    n = len(channels)
    if S.shape != (n, n):
        msg = f"spillover must be {n}x{n} for channels {channels}, got {S.shape}"
        raise ValueError(msg)
    if not np.all(np.isfinite(S)):
        msg = "spillover contains non-finite entries"
        raise ValueError(msg)
    try:
        inv = np.linalg.inv(S)
    except np.linalg.LinAlgError as e:
        msg = "spillover matrix is singular; compensation undefined"
        raise ValueError(msg) from e
    if not np.all(np.isfinite(inv)):
        msg = "spillover inverse produced non-finite values"
        raise ValueError(msg)

    out = events.copy()
    corrected = events[channels].to_numpy(dtype=float) @ inv
    out.loc[:, channels] = corrected
    return out
