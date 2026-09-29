#!/usr/bin/env python3
"""Write deterministic gated-population FCS fixtures into a directory.

Delegates to cytof_qc.demo_fcs (fcsio writer under the hood) so the
pipeline exercises the same fixture path the Snakemake DAG uses.
"""

import importlib.util
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "fcs_io" / "src"))
spec = importlib.util.spec_from_file_location(
    "demo_fcs", REPO / "cytof_qc" / "cytof_qc" / "demo_fcs.py")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

out_dir = sys.argv[1] if len(sys.argv) > 1 else "demo_out"
n_events = int(sys.argv[2]) if len(sys.argv) > 2 else 200
mod.main(out_dir, n_events=n_events)
