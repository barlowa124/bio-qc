#!/usr/bin/env python3
"""Write a run manifest for a bio-qc workflow run.

Common schema (bio-qc/run-manifest@1) emitted by both the Snakemake
and Nextflow paths: engine, git state, params, input/output checksums,
and tool versions. One file per run, next to the results.

Usage:
    python scripts/run_manifest.py --pipeline fcs_qc \
        --engine 'nextflow@26.04.6' --out results/run_manifest.json \
        --inputs data/demo --artifacts results/demo_qc \
        --param demo_events=200
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
from datetime import datetime, timezone


def _sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _git(repo: str) -> dict:
    def g(*args):
        try:
            return subprocess.run(
                ["git", "-C", repo, *args], capture_output=True,
                text=True, timeout=10).stdout.strip()
        except Exception:
            return ""
    return {"sha": g("rev-parse", "HEAD"),
            "dirty": bool(g("status", "--porcelain"))}


def _files(d: str) -> list[dict]:
    out = []
    for root, _, names in os.walk(d):
        for n in sorted(names):
            if n == "run_manifest.json" or n.startswith("."):
                continue
            p = os.path.join(root, n)
            out.append({"path": os.path.relpath(p, d),
                        "sha256": _sha256(p),
                        "bytes": os.path.getsize(p)})
    return sorted(out, key=lambda x: x["path"])


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pipeline", required=True)
    ap.add_argument("--engine", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--inputs")
    ap.add_argument("--artifacts")
    ap.add_argument("--param", action="append", default=[])
    ap.add_argument("--status", default="success")
    a = ap.parse_args()

    repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    manifest = {
        "schema": "bio-qc/run-manifest@1",
        "pipeline": a.pipeline,
        "engine": a.engine,
        "git": _git(repo),
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "status": a.status,
        "params": dict(p.split("=", 1) for p in a.param),
        "inputs": _files(a.inputs) if a.inputs else [],
        "artifacts": _files(a.artifacts) if a.artifacts else [],
        "tools": {"python": sys.version.split()[0]},
    }
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w") as f:
        json.dump(manifest, f, indent=1)
    print(f"wrote {a.out} ({len(manifest['inputs'])} inputs, "
          f"{len(manifest['artifacts'])} artifacts)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
