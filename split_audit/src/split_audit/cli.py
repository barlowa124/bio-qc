"""split-audit CLI: audit a train/test split for leakage.

    split-audit demo                          run the planted-leakage demo
    split-audit scaffolds DATA.csv            exact scaffold-id audit
    split-audit sequences DATA.csv            k-mer Jaccard audit

CSV inputs need an id column, a split column, and a scaffold_id or
sequence column (override with --id-col/--split-col/--scaffold-col/
--seq-col). Writes a JSON result plus a markdown report under --out.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from . import audit, demo, report


def _rows(path: str) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def _emit(result: dict, name: str, kind: str, outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / f"{name}.json").write_text(json.dumps(result, indent=2) + "\n")
    (outdir / f"{name}.md").write_text(report.render(name, result, kind))
    print(report.render(name, result, kind))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="split-audit", description=__doc__)
    sub = parser.add_subparsers(dest="cmd", required=True)

    for cmd in ("scaffolds", "sequences"):
        p = sub.add_parser(cmd)
        p.add_argument("csv")
        p.add_argument("--id-col", default="id")
        p.add_argument("--split-col", default="split")
        p.add_argument("--out", type=Path, default=Path("results"))
    sub.choices["scaffolds"].add_argument("--scaffold-col",
                                          default="scaffold_id")
    seq = sub.choices["sequences"]
    seq.add_argument("--seq-col", default="sequence")
    seq.add_argument("-k", type=int, default=3)
    seq.add_argument("--threshold", type=float, default=0.6)

    p_demo = sub.add_parser("demo")
    p_demo.add_argument("--out", type=Path, default=Path("results"))

    args = parser.parse_args(argv)

    if args.cmd == "demo":
        scaf_rows, seq_rows, planted = demo.make_demo()
        scaf = audit.audit_scaffolds(scaf_rows)
        seqs = audit.audit_sequences(seq_rows)
        _emit(scaf, "scaffolds", "scaffolds", args.out)
        _emit(seqs, "sequences", "sequences", args.out)
        found = set(seqs["flagged_items"]) & set(planted)
        print(f"planted sequence leaks: {len(planted)}, "
              f"recovered: {len(found)}")
        return 0 if found == set(planted) else 1

    rows = _rows(args.csv)
    if args.cmd == "scaffolds":
        result = audit.audit_scaffolds(rows, args.id_col, args.scaffold_col,
                                       args.split_col)
        kind = "scaffolds"
    else:
        result = audit.audit_sequences(rows, args.id_col, args.seq_col,
                                       args.split_col, args.k, args.threshold)
        kind = "sequences"
    _emit(result, Path(args.csv).stem, kind, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
