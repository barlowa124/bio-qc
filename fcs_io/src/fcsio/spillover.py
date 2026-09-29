"""$SPILLOVER extraction: the compensation matrix travels inside the
file's TEXT segment as a flattened keyword — "n,ch1..chn,v11..vnn".
"""

from __future__ import annotations


def spillover(meta: dict) -> tuple[list[str], list[list[float]]] | None:
    """Return (channel_names, matrix) from $SPILLOVER/SPILLOVER, or None."""
    raw = meta.get("$SPILLOVER") or meta.get("SPILLOVER")
    if not raw:
        return None
    parts = [p.strip() for p in raw.split(",")]
    n = int(parts[0])
    names = parts[1:1 + n]
    vals = [float(x) for x in parts[1 + n:1 + n + n * n]]
    if len(vals) != n * n:
        return None
    matrix = [vals[i * n:(i + 1) * n] for i in range(n)]
    return names, matrix
