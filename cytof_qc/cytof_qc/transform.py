"""Signal transforms for mass-cytometry intensity values."""

from __future__ import annotations

import numpy as np
import pandas as pd

DEFAULT_COFACTOR = 5.0


def arcsinh(df: pd.DataFrame, cofactor: float = DEFAULT_COFACTOR) -> pd.DataFrame:
    """Apply the standard CyTOF arcsinh(x / cofactor) transform columnwise.

    Negative inputs (background-subtracted events) transform to negative
    outputs rather than NaN, unlike a log transform.
    """
    if cofactor <= 0:
        msg = f"cofactor must be positive, got {cofactor}"
        raise ValueError(msg)
    return pd.DataFrame(
        np.arcsinh(df.to_numpy(dtype=float) / cofactor),
        index=df.index,
        columns=df.columns,
    )


def inverse_arcsinh(df: pd.DataFrame, cofactor: float = DEFAULT_COFACTOR) -> pd.DataFrame:
    """Invert :func:`arcsinh`; used by tests and for back-translated axes."""
    return pd.DataFrame(
        np.sinh(df.to_numpy(dtype=float)) * cofactor,
        index=df.index,
        columns=df.columns,
    )
