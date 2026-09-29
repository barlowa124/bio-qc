"""FCS loading for mass-cytometry files.

Each FCS file from the gated benchmark sets holds one manually gated
population; the population label comes from the filename convention
``<sample>_<population>_cells.fcs``.
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd
from fcsparser import parse


def load_fcs(path: str | Path) -> tuple[pd.DataFrame, dict]:
    """Read an FCS file and return events renamed to marker names.

    The FCS channel column is named by isotope mass (``$PnS``, e.g. ``115``)
    while the marker name lives in ``$PnN`` (e.g. ``CD45``). Returns the
    event table keyed by marker name, plus the raw metadata dict.
    """
    meta, df = parse(str(path), reformat_meta=False)
    rename = {}
    for k, v in meta.items():
        m = re.match(r"^\$P(\d+)S$", str(k))
        if m:
            channel = str(v)
            marker = meta.get(f"$P{m.group(1)}N")
            if channel in df.columns and marker:
                rename[channel] = str(marker)
    df = df.rename(columns=rename)
    return df, meta


_RESIDUAL_LABELS = {"NotGated", "NotDebrisSinglets", "ungated", "unknown"}


def _parse_filename(name: str) -> tuple[str, str | None]:
    """Parse (population, batch) from a gated-benchmark FCS filename.

    Two conventions are supported:
    - ``Marrow1_Naive CD8+ T_cells.fcs`` (Levine_13dim) -> ("Naive CD8+ T", None)
    - ``..._normalized_Basophils_H1.fcs`` (Levine_32dim) -> ("Basophils", "H1")

    Residual/ungated files normalize to ``NotGated`` so downstream code has
    one label for the ungated remainder regardless of dataset vocabulary.
    """
    m = re.match(r"^.+?_(.+)_cells\.fcs$", name)
    if m:
        pop, batch = m.group(1), None
    else:
        # population is the final underscore-free segment before the
        # batch suffix: "..._normalized_Basophils_H1.fcs" -> "Basophils"
        m = re.match(r"^.+_([^_]+)_(H\d+)\.fcs$", name)
        if not m:
            return "unknown", None
        pop, batch = m.group(1), m.group(2)
    if pop in _RESIDUAL_LABELS:
        pop = "NotGated"
    return pop, batch


def load_population_dir(directory: str | Path) -> pd.DataFrame:
    """Load every gated FCS file under `directory` into one table.

    Adds ``population``, ``batch`` (when the filename carries a donor/run
    suffix), and ``source_file`` columns. Residual ungated files are
    labeled ``NotGated`` across dataset conventions.
    """
    frames = []
    for fcs in sorted(Path(directory).glob("*.fcs")):
        df, _ = load_fcs(fcs)
        pop, batch = _parse_filename(fcs.name)
        df["population"] = pop
        df["batch"] = batch
        df["source_file"] = fcs.name
        frames.append(df)
    if not frames:
        msg = f"no .fcs files under {directory}"
        raise FileNotFoundError(msg)
    return pd.concat(frames, ignore_index=True)
