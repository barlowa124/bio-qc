#!/usr/bin/env python3
"""FCS -> events TSV + metadata JSON, via fcsio (bio-qc/fcs_io).

nf-core module script: PYTHONPATH is set by the process script block so
the sibling package resolves without installation.
"""

import argparse
import json
import sys


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("fcs")
    ap.add_argument("-o", "--events", required=True)
    ap.add_argument("--meta", required=True)
    a = ap.parse_args()

    import fcsio
    events, meta = fcsio.read(a.fcs)
    channels = meta["channels"]
    with open(a.events, "w") as f:
        f.write("\t".join(channels) + "\n")
        for row in events:
            f.write("\t".join(f"{v:.6g}" for v in row) + "\n")
    with open(a.meta, "w") as f:
        json.dump({k: v for k, v in meta.items()
                   if isinstance(v, (str, int, float, list))},
                  f, indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
