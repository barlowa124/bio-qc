"""Event-level and channel-level QC for mass-cytometry event tables.

Flags data-quality problems rather than dropping rows silently: callers
get a report dict they can inspect or gate on.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

DEFAULTS = {
    # fraction of events below zero tolerated before flagging a channel
    "max_negative_fraction": 0.30,
    # fraction of exactly-zero events tolerated (dead detector / unsaturated)
    "max_zero_fraction": 0.60,
    # per-population minimum event count before flagging "sparse"
    "min_population_events": 20,
}


def channel_qc(
    events: pd.DataFrame,
    channels: list[str],
    *,
    max_negative_fraction: float = DEFAULTS["max_negative_fraction"],
    max_zero_fraction: float = DEFAULTS["max_zero_fraction"],
) -> dict[str, dict]:
    """Per-channel rates of negative and zero events, with flag booleans."""
    report = {}
    for ch in channels:
        x = events[ch].to_numpy(dtype=float)
        n = len(x)
        neg = float(np.mean(x < 0)) if n else 0.0
        zero = float(np.mean(x == 0)) if n else 1.0
        report[ch] = {
            "n_events": int(n),
            "negative_fraction": round(neg, 4),
            "zero_fraction": round(zero, 4),
            "median": float(np.median(x)) if n else float("nan"),
            "flag_negative": neg > max_negative_fraction,
            "flag_zero": zero > max_zero_fraction,
        }
    return report


def population_qc(
    events: pd.DataFrame,
    *,
    min_events: int = DEFAULTS["min_population_events"],
    population_col: str = "population",
) -> dict[str, dict]:
    """Per-population event counts, flagging populations too sparse to gate."""
    counts = events[population_col].value_counts()
    return {
        pop: {"n_events": int(n), "flag_sparse": int(n) < min_events}
        for pop, n in counts.items()
    }
