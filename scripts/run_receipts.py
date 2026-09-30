#!/usr/bin/env python3
"""Hash-chained execution receipts for a Nextflow run.

One JSONL per run under results/. Each record binds a task to what it
actually did: process name and tag, the task's staged inputs and written
outputs (sha256), an environment fingerprint (declared versions from the
task's versions.yml plus the engine), the task's own Nextflow hash, and
start/complete timestamps. Records chain via chain_prev like the
inference_receipts format, so a deleted or reordered task breaks the
chain.

Inputs are the task workdir's symlinks (Nextflow stages inputs as
links); outputs are the real files the task wrote, minus Nextflow's
dotfiles. The receipts are a durable record; re-hashing the referenced
files is possible only while work/ or the published outputs survive —
documented under Limitations below and in the README.

Usage:
    python scripts/run_receipts.py --trace results/trace.txt \
        --out results/run_receipts.jsonl --run-name sad_euler \
        --session-id 4f9a... --engine nextflow@26.04.6

    python scripts/run_receipts.py --verify results/run_receipts.jsonl
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import sys
import time
from datetime import datetime, timezone

SCHEMA_VERSION = 1

# Nextflow bookkeeping files inside a task workdir — not outputs.
_NF_DOTFILES = (".command.", ".exitcode", ".nf-", ".nextflow")


def canonical_json(obj) -> bytes:
    return json.dumps(obj, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False).encode("utf-8")


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _hash_dir_entries(workdir: str) -> tuple[list[dict], list[dict]]:
    """Inputs = symlinks (staged); outputs = real files written."""
    inputs, outputs = [], []
    if not os.path.isdir(workdir):
        return inputs, outputs
    for name in sorted(os.listdir(workdir)):
        if name.startswith("."):
            continue  # .command.*, .exitcode, .begin etc.
        p = os.path.join(workdir, name)
        if os.path.islink(p):
            target = os.path.realpath(p)
            inputs.append({"name": name, "staged_from": target,
                           "sha256": sha256_file(target)})
        elif os.path.isfile(p):
            outputs.append({"name": name, "sha256": sha256_file(p)})
    return inputs, outputs


def _env_fingerprint(workdir: str, process: str, engine: str) -> dict:
    """Declared tool versions (versions.yml) + engine identity."""
    env = {"process": process, "engine": engine}
    vyml = os.path.join(workdir, "versions.yml")
    if os.path.isfile(vyml):
        env["versions_yml"] = open(vyml, encoding="utf-8").read()
    return env


def _env_hash(env: dict) -> str:
    return hashlib.sha256(canonical_json(env)).hexdigest()


def _stable_trace(path: str, tries: int = 25, interval: float = 0.2) -> None:
    """Wait for the trace file to stop growing — the writer may not have
    flushed the last row when the workflow tail runs."""
    last = -1
    for _ in range(tries):
        try:
            size = os.path.getsize(path)
        except FileNotFoundError:
            size = -1
        if size == last and size > 0:
            return
        last = size
        time.sleep(interval)


def load_trace(path: str) -> list[dict]:
    with open(path, newline="", encoding="utf-8") as f:
        return [r for r in csv.DictReader(f, delimiter="\t")]


def task_receipts(trace_path: str, *, engine: str, run_name: str,
                  session_id: str) -> list[dict]:
    rows = load_trace(trace_path)
    rows.sort(key=lambda r: r.get("submit") or r.get("start") or "")
    out, prev_id = [], "genesis"
    for r in rows:
        workdir = r.get("workdir", "")
        process = r.get("process") or r.get("name", "").split(" ")[0]
        inputs, outputs = _hash_dir_entries(workdir)
        env = _env_fingerprint(workdir, process, engine)
        body = {
            "v": SCHEMA_VERSION,
            "kind": "task_receipt",
            "run_name": run_name,
            "session_id": session_id,
            "task": r.get("name", ""),
            "process": process,
            "workdir": workdir,
            "status": r.get("status", ""),
            "exit": r.get("exit", ""),
            "nextflow_task_hash": r.get("hash", ""),
            "env": env,
            "env_sha256": _env_hash(env),
            "inputs": inputs,
            "outputs": outputs,
            "started_at": r.get("start", ""),
            "completed_at": r.get("complete", ""),
            "duration_s": r.get("realtime", ""),
            "recorded_at": datetime.now(timezone.utc).isoformat(),
            "chain_prev": prev_id,
        }
        rid = "rcpt-" + hashlib.sha256(canonical_json(body)).hexdigest()[:16]
        body["receipt_id"] = rid
        prev_id = rid
        out.append(body)
    return out


def _receipt_hash(r: dict) -> str:
    body = {k: v for k, v in r.items() if k != "receipt_id"}
    return "rcpt-" + hashlib.sha256(canonical_json(body)).hexdigest()[:16]


def check_chain(receipts: list[dict]) -> list[str]:
    """Integrity problems in a receipt log (tampered body, chain break)."""
    problems = []
    prev_id = "genesis"
    for i, r in enumerate(receipts):
        if _receipt_hash(r) != r.get("receipt_id"):
            problems.append(f"[{i}] {r.get('receipt_id', '?')}: content "
                            "hash does not match receipt_id")
        if r.get("env_sha256") != _env_hash(r.get("env", {})):
            problems.append(f"[{i}] {r.get('receipt_id', '?')}: env hash "
                            "does not match env record")
        if r.get("chain_prev") != prev_id:
            problems.append(f"[{i}] {r.get('receipt_id', '?')}: chain "
                            f"break (expected prev={prev_id})")
        prev_id = r.get("receipt_id", "?")
    return problems


def rehash(receipts: list[dict]) -> list[str]:
    """Re-hash the recorded files where they still exist.

    Chain verification proves the record is intact; this proves the
    files still match the record. Inputs resolve through the recorded
    staged symlink target, outputs through workdir/name. Files that no
    longer exist are reported as unavailable, not as mismatches — the
    distinction matters when work/ has been cleaned.
    """
    problems, checked = [], 0
    for i, r in enumerate(receipts):
        for e in r.get("inputs", []):
            p = e.get("staged_from", "")
            if p and os.path.exists(p):
                checked += 1
                if sha256_file(p) != e.get("sha256"):
                    problems.append(f"[{i}] {r.get('task')}: input "
                                    f"{e.get('name')} hash mismatch")
        wd = r.get("workdir", "")
        for e in r.get("outputs", []):
            p = os.path.join(wd, e.get("name", "")) if wd else ""
            if p and os.path.exists(p):
                checked += 1
                if sha256_file(p) != e.get("sha256"):
                    problems.append(f"[{i}] {r.get('task')}: output "
                                    f"{e.get('name')} hash mismatch")
    if not checked:
        problems.append("no recorded files still exist — record "
                        "integrity is all that can be verified")
    else:
        problems.insert(0, f"re-hashed {checked} recorded files")
    return problems


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--trace", help="nextflow -with-trace TSV")
    ap.add_argument("--out", help="receipts JSONL path")
    ap.add_argument("--run-name", default="")
    ap.add_argument("--session-id", default="")
    ap.add_argument("--engine", default="nextflow")
    ap.add_argument("--verify", help="check an existing receipts file")
    ap.add_argument("--rehash", help="re-hash recorded files that "
                    "still exist and compare to the receipts")
    a = ap.parse_args()

    if a.verify:
        receipts = [json.loads(l) for l in open(a.verify) if l.strip()]
        problems = check_chain(receipts)
        print(f"{a.verify}: {len(receipts)} receipts, "
              f"{len(problems)} problems")
        for p in problems:
            print("  " + p)
        return 1 if problems else 0

    if a.rehash:
        receipts = [json.loads(l) for l in open(a.rehash) if l.strip()]
        problems = check_chain(receipts)
        print(f"{a.rehash}: {len(receipts)} receipts, "
              f"{len(problems)} chain problems")
        for p in problems:
            print("  " + p)
        issues = rehash(receipts)
        for p in issues:
            print("  " + p)
        return 1 if problems or any("mismatch" in i for i in issues) \
            else 0

    if not (a.trace and a.out):
        ap.error("--trace and --out required (or use --verify)")
    _stable_trace(a.trace)
    receipts = task_receipts(a.trace, engine=a.engine,
                             run_name=a.run_name,
                             session_id=a.session_id)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w", encoding="utf-8") as f:
        for r in receipts:
            f.write(json.dumps(r, ensure_ascii=False) + "\n")
    print(f"receipts: {len(receipts)} task records -> {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
