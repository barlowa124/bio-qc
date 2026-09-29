"""FCS 3.x reader. See package docstring for the quirk policy.

Layout per spec: a 58-byte HEADER (version + six 8-char ASCII offset
fields), a TEXT segment of delimiter-separated key/value pairs, then the
binary DATA segment described by TEXT keywords.
"""

from __future__ import annotations

import struct


class FCSError(ValueError):
    """Malformed or unsupported FCS file."""


def _ascii_offsets(header: bytes) -> tuple[int, int, int, int]:
    """Six 8-char right-justified fields: TEXT begin/end, DATA begin/end,
    ANALYSIS begin/end (unused). Blank or '0' fields mean the offset lives
    in a TEXT keyword instead."""
    fields = []
    for i in range(6):
        raw = header[10 + i * 8: 18 + i * 8].strip()
        fields.append(int(raw) if raw else 0)
    return fields[0], fields[1], fields[2], fields[3]


def _read_text(raw: bytes, t0: int, t1: int) -> dict:
    """Split the TEXT segment on its delimiter; doubled delimiters inside
    a value are escaped literals and get merged back via a sentinel."""
    seg = raw[t0:t1 + 1]
    if not seg:
        raise FCSError("empty TEXT segment")
    delim = seg[0]
    body = seg[1:]
    SENT = b"\x00"
    body = body.replace(bytes([delim, delim]), SENT)
    parts = body.split(bytes([delim]))
    meta: dict[str, str] = {}
    it = iter(parts)
    for k in it:
        v = next(it, b"")
        key = k.decode("utf-8", "replace").strip()
        val = v.replace(SENT, bytes([delim])).decode(
            "utf-8", "replace").strip()
        if key:
            meta[key] = val
    return meta


def _byte_order(byteord: str) -> str:
    """Map $BYTEORD ("1,2,3,4" = least->most significant) to struct."""
    order = [int(x) for x in byteord.split(",")]
    if order[:2] == [1, 2]:
        return "<"          # little-endian families (1,2,3,4 / 1,2)
    if order[:2] == [4, 3] or order == sorted(order, reverse=True):
        return ">"          # big-endian (4,3,2,1 / 2,1)
    raise FCSError(f"unsupported $BYTEORD {byteord!r}")


def _read_data(raw: bytes, d0: int, d1: int, meta: dict) -> list[list]:
    datatype = meta.get("$DATATYPE", "F").upper()
    n_par = int(meta.get("$PAR", "0"))
    n_tot = int(meta.get("$TOT", "0"))
    if n_par <= 0 or n_tot < 0:
        raise FCSError("$PAR/$TOT missing or invalid")
    seg = raw[d0:d1 + 1]
    bo = _byte_order(meta.get("$BYTEORD", "1,2,3,4"))

    if datatype == "I":
        widths = {int(meta[f"$P{i + 1}B"]) for i in range(n_par)}
        if len(widths) != 1:
            raise FCSError(
                "mixed $PnB bit widths unsupported "
                f"({sorted(widths)}) — vendor quirk, not guessed")
        width = widths.pop()
        fmt_char = {8: "B", 16: "H", 32: "I"}.get(width)
        if fmt_char is None:
            raise FCSError(f"unsupported integer width {width}")
        row = struct.Struct(bo + fmt_char * n_par)
    elif datatype in ("F", "D"):
        fmt_char = "f" if datatype == "F" else "d"
        row = struct.Struct(bo + fmt_char * n_par)
    else:
        raise FCSError(f"unsupported $DATATYPE {datatype!r}")

    needed = row.size * n_tot
    if len(seg) < needed:
        raise FCSError(f"DATA short: need {needed}, have {len(seg)}")
    return [list(row.unpack_from(seg, i * row.size))
            for i in range(n_tot)]


def read(path: str) -> tuple[list[list], dict]:
    """Read an FCS file. Returns (events, meta): events is a list of
    n_tot rows of n_par values, meta the TEXT keyword dict plus derived
    fields (channels, stains, spillover in spillover.py, nextdata)."""
    raw = open(path, "rb").read()
    if len(raw) < 58 or raw[:3] != b"FCS":
        raise FCSError("not an FCS file (missing FCSx.y magic)")
    meta = _read_text(raw, *_ascii_offsets(raw)[:2])
    meta["$VERSION"] = raw[:6].decode("ascii", "replace")

    _, _, d0, d1 = _ascii_offsets(raw)
    if d0 == 0:
        # Header left the offsets blank — the values live in keywords.
        d0 = int(meta.get("$BEGINDATA", "0"))
        # $ENDDATA is the offset of the last data byte (inclusive), same
        # convention as the header field — do not adjust it.
        d1 = int(meta.get("$ENDDATA", "0"))
    events = _read_data(raw, d0, d1, meta)

    n_par = int(meta["$PAR"])
    meta["channels"] = [meta.get(f"$P{i + 1}N", f"P{i + 1}")
                        for i in range(n_par)]
    meta["stains"] = [meta.get(f"$P{i + 1}S", "") for i in range(n_par)]
    meta["n_events"] = len(events)
    meta["nextdata"] = int(meta.get("$NEXTDATA", "0") or 0)
    return events, meta
