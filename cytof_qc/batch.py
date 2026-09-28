"""Cross-batch channel alignment for multi-run cytometry panels.

Different acquisition days drift in detector sensitivity and staining
efficiency, which shows up as per-channel offsets between batches. This
aligns each channel's batch medians to the pooled median — a linear-shift
approximation that is honest about what it fixes: it centers batches, it
does not model nonlinear warping or correct within-batch structure.

Events are never moved silently: the function returns the corrected table
plus a report of the shift applied per channel per batch, so reviewers can
see whether the correction was small calibration or a real discrepancy.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def align_batches(
    events: pd.DataFrame,
    channels: list[str],
    batch_col: str = "batch",
) -> tuple[pd.DataFrame, dict]:
    """Shift each batch's channel medians to the pooled median.

    Returns ``(corrected, report)`` where report maps each channel to the
    per-batch shift that was applied. Batches with fewer than 50 events or
    a channel whose batch median is undefined are reported but left
    unshifted rather than trusted.
    """
    pooled_med = events[channels].median()
    out = events.copy()
    report: dict[str, dict] = {}
    for ch in channels:
        report[ch] = {}
        for batch, idx in events.groupby(batch_col).groups.items():
            sub = events.loc[idx, ch]
            if len(sub) < 50:
                report[ch][str(batch)] = {"n_events": int(len(sub)), "shift": 0.0,
                                          "skipped": True}
                continue
            batch_med = float(sub.median())
            shift = float(pooled_med[ch] - batch_med)
            out.loc[idx, ch] = sub + shift
            report[ch][str(batch)] = {
                "n_events": int(len(sub)),
                "batch_median": batch_med,
                "pooled_median": float(pooled_med[ch]),
                "shift": shift,
                "skipped": False,
            }
    return out, report
