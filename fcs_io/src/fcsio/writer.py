"""FCS 3.0 writer — minimal valid files, primarily for fixtures and demos.

Writes float32 little-endian DATA with header offsets; the reader module
round-trips what this emits. TEXT values get delimiter-escaped per spec.
"""

from __future__ import annotations

import struct

DELIM = b"\x0c"


def write(path: str, events: list[list], meta: dict) -> None:
    """Serialize events+metadata as FCS3.0. `meta` supplies extra TEXT
    keywords; $PAR/$TOT/$BYTEORD/$DATATYPE/$PnN/$PnR/$PnE/$MODE are
    derived here and must not be overridden."""
    n_par = len(events[0]) if events else len(meta.get("channels", []))
    body = {
        "$MODE": "L", "$BYTEORD": "1,2,3,4", "$DATATYPE": "F",
        "$TOT": str(len(events)), "$PAR": str(n_par),
        "$NEXTDATA": "0",
    }
    body.update({k: str(v) for k, v in meta.items()
                 if not k.startswith("$P")})
    channels = meta.get("channels") or [f"P{i+1}" for i in range(n_par)]
    stains = meta.get("stains") or [""] * n_par
    for i in range(n_par):
        body[f"$P{i+1}N"] = channels[i]
        body[f"$P{i+1}B"] = "32"
        body[f"$P{i+1}R"] = "262144"
        body[f"$P{i+1}E"] = "0,0"
        if stains[i]:
            body[f"$P{i+1}S"] = stains[i]

    def esc(v: str) -> bytes:
        return v.encode().replace(DELIM, DELIM * 2)

    text = DELIM + DELIM.join(
        esc(k) + DELIM + esc(v) for k, v in body.items()) + DELIM

    row = struct.Struct("<" + "f" * n_par)
    data = b"".join(row.pack(*[float(x) for x in e]) for e in events)

    header_len = 58
    t0, t1 = header_len, header_len + len(text) - 1
    d0, d1 = t1 + 1, t1 + len(data)
    header = (b"FCS3.0".ljust(10)
              + f"{t0:8d}{t1:8d}{d0:8d}{d1:8d}{0:8d}{0:8d}".encode())
    with open(path, "wb") as f:
        f.write(header + text + data)
