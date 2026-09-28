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


def load_population_dir(directory: str | Path) -> pd.DataFrame:
    """Load every ``*_cells.fcs`` file under `directory` into one table.

    Adds a ``population`` column parsed from the filename
    (``Marrow1_Naive CD8+ T_cells.fcs`` -> ``Naive CD8+ T``). Files whose
    names do not match the convention are labeled ``unknown``.
    """
    frames = []
    for fcs in sorted(Path(directory).glob("*.fcs")):
        df, _ = load_fcs(fcs)
        m = re.match(r"^.+?_(.+)_cells\.fcs$", fcs.name)
        df["population"] = m.group(1) if m else "unknown"
        df["source_file"] = fcs.name
        frames.append(df)
    if not frames:
        msg = f"no .fcs files under {directory}"
        raise FileNotFoundError(msg)
    return pd.concat(frames, ignore_index=True)
