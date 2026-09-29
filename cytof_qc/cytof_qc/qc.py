"""Event-level and channel-level QC for mass-cytometry event tables.

Flags data-quality problems rather than dropping rows silently: callers
get a report dict they can inspect or gate on.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

# Channels that instrument the run rather than phenotype the cell. They
# are reported by QC but must never drive clustering: DNA channels mark
# cell events (and carry bead-normalization drift), viability is a
# live/dead gate, and Time/event_number are acquisition bookkeeping.
QC_ONLY_CHANNEL_NAMES = {"Time", "event_number", "Event_length", "beadDist"}
QC_ONLY_CHANNEL_PREFIXES = ("DNA", "Viability", "Bead")


def split_channels(channels: list[str]) -> tuple[list[str], list[str]]:
    """Partition columns into phenotypic markers and QC-only channels."""
    pheno, qc_only = [], []
    for ch in channels:
        if ch in QC_ONLY_CHANNEL_NAMES or ch.startswith(
            QC_ONLY_CHANNEL_PREFIXES
        ):
            qc_only.append(ch)
        else:
            pheno.append(ch)
    return pheno, qc_only


DEFAULTS = {
    # fraction of events below zero tolerated before flagging a channel
    "max_negative_fraction": 0.30,
    # fraction of exactly-zero events tolerated (dead detector / unsaturated)
    "max_zero_fraction": 0.60,
    # per-population minimum event count before flagging "sparse"
    "min_population_events": 20,
    # acquisition-drift gate: bins of event order, and the allowed spread
    # of per-bin channel medians in arcsinh units before flagging drift.
    # 0.5 arcsinh units at cofactor 5 is roughly a 1.6x intensity change.
    "drift_bins": 10,
    "max_median_drift": 0.5,
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


def acquisition_qc(
    events: pd.DataFrame,
    channels: list[str],
    *,
    n_bins: int = DEFAULTS["drift_bins"],
    max_median_drift: float = DEFAULTS["max_median_drift"],
    batch_col: str | None = None,
) -> dict:
    """Flag channels whose signal drifts over acquisition order.

    Events are binned by row order (a proxy for acquisition time on a
    sequential instrument) and each channel's per-bin median is tracked.
    A channel is flagged when the spread of its bin medians exceeds
    ``max_median_drift`` — the signature of bead-normalization failure or
    detector decay mid-run.

    When ``batch_col`` is set, drift is measured within each batch
    separately and the report is keyed ``{batch: {channel: ...}}``;
    otherwise the batch seam would masquerade as in-run drift.
    """
    if batch_col is not None and batch_col in events.columns:
        return {
            str(b): acquisition_qc(
                sub.reset_index(drop=True), channels,
                n_bins=n_bins, max_median_drift=max_median_drift,
            )
            for b, sub in events.groupby(batch_col)
        }
    report = {}
    n = len(events)
    bins = np.array_split(np.arange(n), n_bins) if n else []
    for ch in channels:
        x = events[ch].to_numpy(dtype=float)
        meds = np.array([np.median(x[b]) for b in bins], dtype=float) if bins else np.array([np.nan])
        drift = float(np.nanmax(meds) - np.nanmin(meds)) if len(meds) else float("nan")
        report[ch] = {
            "bin_medians": [round(float(v), 4) for v in meds],
            "median_drift": round(drift, 4),
            "flag_drift": bool(np.isfinite(drift) and drift > max_median_drift),
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
