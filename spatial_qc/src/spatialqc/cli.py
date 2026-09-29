"""spatialqc — Visium QC metrics and filter-strategy comparison.

  python -m spatialqc.cli SAMPLE_DIR --out results/spatial_qc.json
                       [--no-cluster]

SAMPLE_DIR points at a Space Ranger export (filtered_feature_bc_matrix/
+ spatial/tissue_positions_list.csv), or at the matrix dir itself.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="spatialqc", description=__doc__)
    ap.add_argument("sample_dir")
    ap.add_argument("--out", default=None)
    ap.add_argument("--dataset", default=None,
                    help="dataset label recorded in the report")
    ap.add_argument("--no-cluster", action="store_true")
    args = ap.parse_args(argv)

    from . import io, metrics
    from .compare import compare as run_compare
    sample = io.load_visium(args.sample_dir)
    met = metrics.spot_metrics(sample)
    report = {
        "dataset": args.dataset or Path(args.sample_dir).name,
        "summary": metrics.summary(sample, met),
        "strategy_comparison": run_compare(
            sample, run_clusters=not args.no_cluster),
    }
    text = json.dumps(report, indent=1)
    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(text + "\n")
        s = report["summary"]
        print(f"{s['n_spots']} spots ({s['n_in_tissue']} in tissue), "
              f"{s['n_genes']} genes -> {args.out}", file=sys.stderr)
    else:
        print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
