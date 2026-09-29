#!/usr/bin/env python3
"""Summarize service_requests.jsonl — the read half of the instrument
loop. Prints request counts, error rate, latency percentiles, and upload
sizes so the next change responds to measured behavior, not guesses."""

import json
import sys
from collections import Counter
from pathlib import Path

LOG = Path(__file__).resolve().parent.parent / "service_requests.jsonl"


def pctile(sorted_vals, q):
    if not sorted_vals:
        return 0.0
    idx = min(len(sorted_vals) - 1, int(q * len(sorted_vals)))
    return sorted_vals[idx]


def main(path=LOG):
    entries = []
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if line:
            entries.append(json.loads(line))
    if not entries:
        print("no requests logged")
        return

    by_path = Counter(e["path"] for e in entries)
    errors = [e for e in entries if e["status"] >= 400]
    lat = sorted(e["ms"] for e in entries)
    sizes = sorted(e["n_events"] for e in entries if "n_events" in e)

    print(f"{len(entries)} requests, {len(errors)} errors "
          f"({100 * len(errors) / len(entries):.0f}%)")
    for p, n in by_path.most_common():
        print(f"  {p}: {n}")
    print(f"latency ms: p50={pctile(lat, 0.5)} p95={pctile(lat, 0.95)} "
          f"max={lat[-1]}")
    if sizes:
        print(f"upload events: p50={pctile(sizes, 0.5)} "
              f"p95={pctile(sizes, 0.95)} max={sizes[-1]}")
    if errors:
        print("error breakdown:", dict(Counter(e['status'] for e in errors)))


if __name__ == "__main__":
    main(Path(sys.argv[1]) if len(sys.argv) > 1 else LOG)
