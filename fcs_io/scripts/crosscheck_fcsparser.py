"""Cross-validate fcsio against fcsparser on the Levine benchmark FCS set.

Reads every ``*.fcs`` under the given directory twice — once with fcsio,
once with fcsparser — and reports per-file event-matrix agreement plus
channel-name ($PnN) agreement. Writes a compact JSON report; raw FCS
data itself stays out of git.

Usage:
    python scripts/crosscheck_fcsparser.py <fcs_dir> [more_dirs] \
        --out validation/crosscheck_levine.json

Exits nonzero if any file mismatches beyond float tolerance.
"""

from __future__ import annotations

import glob
import json
import sys
from pathlib import Path


def run(dirs: list[str]) -> dict:
    import numpy as np
    from fcsparser import parse as ref_parse
    import fcsio

    files = sorted(
        p for d in dirs for p in glob.glob(str(Path(d) / "*.fcs")))
    report = {"n_files": len(files), "n_events": 0,
              "max_abs_diff": 0.0, "n_exact": 0, "n_allclose": 0,
              "channels_ok": 0, "mismatches": []}
    for path in files:
        ev, meta = fcsio.read(path)
        mref, dref = ref_parse(path, reformat_meta=False)
        mine = np.asarray(ev, dtype=float)
        ref = dref.to_numpy(dtype=float)
        report["n_events"] += len(ev)
        if mine.shape != ref.shape:
            report["mismatches"].append(
                {Path(path).name:
                 f"shape {mine.shape} vs ref {ref.shape}"})
            continue
        diff = float(np.abs(mine - ref).max())
        report["max_abs_diff"] = max(report["max_abs_diff"], diff)
        if diff == 0:
            report["n_exact"] += 1
        elif np.allclose(mine, ref):
            report["n_allclose"] += 1
        else:
            report["mismatches"].append(
                {Path(path).name: f"max_abs_diff {diff:.4g}"})
        refn = [str(mref.get(f"$P{i + 1}N", ""))
                for i in range(len(dref.columns))]
        if refn == meta["channels"]:
            report["channels_ok"] += 1
    report["ok"] = (not report["mismatches"]
                    and report["channels_ok"] == report["n_files"])
    return report


if __name__ == "__main__":
    argv = sys.argv[1:]
    out = argv[argv.index("--out") + 1] if "--out" in argv else None
    dirs = [a for a in argv if not a.startswith("-") and a != out]
    rep = run(dirs)
    text = json.dumps(rep, indent=1)
    if out:
        Path(out).parent.mkdir(parents=True, exist_ok=True)
        Path(out).write_text(text + "\n")
        print(f"wrote {out}")
    else:
        print(text)
    sys.exit(0 if rep["ok"] else 1)
