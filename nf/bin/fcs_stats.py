#!/usr/bin/env python3
"""Per-channel summary stats over an events TSV -> JSON.

Counts, per-channel mean/median/p95, negative fraction — the QC
surface a flow bench looks at first. Stdlib only.
"""

import argparse
import json
import statistics
import sys


def _pct(vals: list[float], q: float) -> float:
    s = sorted(vals)
    k = (len(s) - 1) * (q / 100)
    lo, hi = int(k), min(int(k) + 1, len(s) - 1)
    return s[lo] + (s[hi] - s[lo]) * (k - lo)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("events_tsv")
    ap.add_argument("-o", "--out", required=True)
    a = ap.parse_args()

    with open(a.events_tsv) as f:
        header = f.readline().rstrip("\n").split("\t")
        cols = {i: [] for i in range(len(header))}
        n = 0
        for line in f:
            n += 1
            for i, v in enumerate(line.rstrip("\n").split("\t")):
                cols[i].append(float(v))
    stats = {}
    for i, name in enumerate(header):
        v = cols[i]
        stats[name] = {
            "mean": round(statistics.fmean(v), 4),
            "median": round(statistics.median(v), 4),
            "p95": round(_pct(v, 95), 4),
            "neg_frac": round(sum(1 for x in v if x < 0) / n, 4),
        }
    with open(a.out, "w") as f:
        json.dump({"n_events": n, "channels": stats}, f, indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
