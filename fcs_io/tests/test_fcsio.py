"""fcsio tests: build minimal FCS files in-memory and parse them back.

The fixture writer is itself evidence the format is understood — each
test exercises one vendor quirk documented in the package docstring.
"""
import os
import struct
import tempfile
import unittest

from fcsio import read, FCSError, spillover

DELIM = b"\x0c"


def make_fcs(meta_pairs: dict, events: list[list], *, datatype="F",
             byteord="1,2,3,4", version="FCS3.0",
             header_offsets=True) -> bytes:
    """Serialize a minimal valid FCS file."""
    body_meta = dict(meta_pairs)
    n_par = len(events[0]) if events else int(meta_pairs.get("$PAR", 2))
    n_tot = len(events)
    body_meta.setdefault("$PAR", str(n_par))
    body_meta.setdefault("$TOT", str(n_tot))
    body_meta["$BYTEORD"] = byteord
    body_meta["$DATATYPE"] = datatype
    body_meta.setdefault("$MODE", "L")
    for i in range(n_par):
        body_meta.setdefault(f"$P{i+1}N", f"CH{i+1}")
        if datatype == "I":
            body_meta.setdefault(f"$P{i+1}B", "16")
        body_meta.setdefault(f"$P{i+1}R", "262144")
        body_meta.setdefault(f"$P{i+1}E", "0,0")

    bo = "<" if byteord.split(",")[0] == "1" else ">"
    if datatype == "I":
        row = struct.Struct(bo + "H" * n_par)
    else:
        row = struct.Struct(bo + ("f" if datatype == "F" else "d") * n_par)
    data = b"".join(row.pack(*e) for e in events)

    def esc(v: str) -> bytes:
        # spec escape: a delimiter inside a key/value is written doubled
        return v.encode().replace(DELIM, DELIM * 2)

    text = DELIM + DELIM.join(
        esc(k) + DELIM + esc(v)
        for k, v in body_meta.items()) + DELIM

    header_len = 58
    if not header_offsets:
        # blank DATA offsets in header; values ride in $BEGINDATA/
        # $ENDDATA keywords. The keyword text contains the offsets, so
        # solve to a fixed point (digit widths can shift the answer).
        d0, d1 = 0, 0
        while True:
            extra = DELIM.join(
                x.encode() for x in
                ("$BEGINDATA", str(d0), "$ENDDATA", str(d1)))
            t = text + extra + DELIM
            t1 = header_len + len(t) - 1
            nd0, nd1 = t1 + 1, t1 + len(data)
            if (nd0, nd1) == (d0, d1):
                text = t
                break
            d0, d1 = nd0, nd1
        header = (version.encode().ljust(10)
                  + f"{header_len:8d}{t1:8d}".encode() + b" " * 32)
        return header + text + data
    t0, t1 = header_len, header_len + len(text) - 1
    d0, d1 = t1 + 1, t1 + len(data)
    header = (version.encode().ljust(10)
              + f"{t0:8d}{t1:8d}{d0:8d}{d1:8d}{0:8d}{0:8d}".encode())
    return header + text + data


def write(tmp, blob) -> str:
    path = os.path.join(tmp, "test.fcs")
    with open(path, "wb") as f:
        f.write(blob)
    return path


class ReaderTests(unittest.TestCase):
    def test_float32_little_endian(self):
        with tempfile.TemporaryDirectory() as d:
            blob = make_fcs({}, [[1.5, 2.5], [3.0, 4.25]])
            ev, meta = read(write(d, blob))
        self.assertEqual(len(ev), 2)
        self.assertAlmostEqual(ev[0][0], 1.5)
        self.assertEqual(meta["channels"], ["CH1", "CH2"])
        self.assertEqual(meta["n_events"], 2)

    def test_float32_big_endian(self):
        with tempfile.TemporaryDirectory() as d:
            blob = make_fcs({}, [[1.5, 2.5]], byteord="4,3,2,1")
            ev, _ = read(write(d, blob))
        self.assertAlmostEqual(ev[0][1], 2.5)

    def test_integer_datatype(self):
        with tempfile.TemporaryDirectory() as d:
            blob = make_fcs({}, [[100, 200], [65535, 0]], datatype="I")
            ev, meta = read(write(d, blob))
        self.assertEqual(ev[1][0], 65535)
        self.assertEqual(ev[0][1], 200)

    def test_channel_and_stain_names(self):
        with tempfile.TemporaryDirectory() as d:
            blob = make_fcs({"$P1S": "CD45", "$P2S": "CD3"}, [[1.0, 2.0]])
            _, meta = read(write(d, blob))
        self.assertEqual(meta["stains"], ["CD45", "CD3"])

    def test_escaped_delimiter_in_value(self):
        with tempfile.TemporaryDirectory() as d:
            # vendor quirk: a value containing the delimiter escapes it
            # doubled — "a\x0c\x0cb" stores 'a\x0cb'
            meta_pairs = {"CUSTOM": "left\x0cright"}
            blob = make_fcs(meta_pairs, [[1.0, 2.0]])
            _, meta = read(write(d, blob))
        self.assertEqual(meta["CUSTOM"], "left\x0cright")

    def test_keyword_offsets_when_header_blank(self):
        with tempfile.TemporaryDirectory() as d:
            blob = make_fcs({}, [[9.0, 8.0]], header_offsets=False)
            ev, meta = read(write(d, blob))
        self.assertAlmostEqual(ev[0][0], 9.0)

    def test_spillover_matrix(self):
        with tempfile.TemporaryDirectory() as d:
            blob = make_fcs(
                {"$SPILLOVER": "2,CH1,CH2,1.0,0.15,0.02,1.0"},
                [[1.0, 2.0]])
            _, meta = read(write(d, blob))
        names, m = spillover(meta)
        self.assertEqual(names, ["CH1", "CH2"])
        self.assertAlmostEqual(m[0][1], 0.15)

    def test_not_fcs_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            path = write(d, b"NOTANFCSFILE" + b"\x00" * 100)
            with self.assertRaises(FCSError):
                read(path)

    def test_mixed_bitwidths_refused(self):
        with tempfile.TemporaryDirectory() as d:
            blob = make_fcs({"$P2B": "32"}, [[1, 2]], datatype="I")
            with self.assertRaises(FCSError):
                read(write(d, blob))

    def test_nextdata_surfaced(self):
        with tempfile.TemporaryDirectory() as d:
            blob = make_fcs({"$NEXTDATA": "4096"}, [[1.0, 2.0]])
            _, meta = read(write(d, blob))
        self.assertEqual(meta["nextdata"], 4096)


if __name__ == "__main__":
    unittest.main()


class WriterTests(unittest.TestCase):
    def test_roundtrip(self):
        from fcsio import write as fcs_write
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "rt.fcs")
            fcs_write(path, [[1.5, 2.5], [3.0, 4.25]],
                      {"channels": ["FSC-A", "CD45"],
                       "stains": ["", "CD45"]})
            ev, meta = read(path)
        self.assertEqual(meta["channels"], ["FSC-A", "CD45"])
        self.assertEqual(meta["stains"][1], "CD45")
        self.assertEqual(len(ev), 2)
        self.assertAlmostEqual(ev[1][1], 4.25)

    def test_roundtrip_escaped_value(self):
        from fcsio import write as fcs_write
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "rt2.fcs")
            fcs_write(path, [[1.0]],
                      {"channels": ["A"],
                       "COMMENT": "has\x0cdelimiter"})
            _, meta = read(path)
        self.assertEqual(meta["COMMENT"], "has\x0cdelimiter")
