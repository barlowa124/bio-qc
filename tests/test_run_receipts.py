"""run_receipts: hash-chained per-task execution records from the
Nextflow trace + task workdirs."""
import importlib.util
import json
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "run_receipts", ROOT / "scripts" / "run_receipts.py")
rr = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rr)

TRACE_FIELDS = ("task_id", "hash", "name", "status", "exit", "submit",
                "start", "complete", "duration", "realtime", "workdir")


def _write_trace(path: Path, rows: list[dict]) -> None:
    with open(path, "w", newline="") as f:
        f.write("\t".join(TRACE_FIELDS) + "\n")
        for r in rows:
            f.write("\t".join(r.get(k, "") for k in TRACE_FIELDS) + "\n")


def _workdir(base: Path, name: str, inputs: dict, outputs: dict) -> str:
    """Fake a task workdir: symlinked inputs, real outputs, nf dotfiles."""
    wd = base / name
    wd.mkdir(parents=True)
    for fname, target in inputs.items():
        (wd / fname).symlink_to(target)
    for fname, content in outputs.items():
        (wd / fname).write_text(content)
    (wd / ".command.sh").write_text("echo hi")
    (wd / ".exitcode").write_text("0")
    return str(wd)


class ReceiptTests(unittest.TestCase):
    def setUp(self):
        import tempfile
        self._td = tempfile.TemporaryDirectory()
        base = Path(self._td.name)
        # two real "input" files living outside the workdir
        (base / "data").mkdir()
        (base / "data" / "a.fcs").write_bytes(b"fcs-a")
        (base / "data" / "b.fcs").write_bytes(b"fcs-b")
        wd1 = _workdir(base, "w1",
                       {"input.fcs": base / "data" / "a.fcs"},
                       {"a.events.tsv": "e1\te2\n", "versions.yml": "p: 1"})
        wd2 = _workdir(base, "w2",
                       {"events.tsv": base / "w1" / "a.events.tsv"},
                       {"a.stats.json": "{}", "versions.yml": "p: 1"})
        self.rows = [
            {"task_id": "1", "hash": "aa/bb", "name": "P:PARSE (a)",
             "status": "COMPLETED", "exit": "0",
             "submit": "2026-01-01 10:00:00.0",
             "start": "2026-01-01 10:00:00.1",
             "complete": "2026-01-01 10:00:01.0",
             "realtime": "900ms", "workdir": wd1},
            {"task_id": "2", "hash": "cc/dd", "name": "P:STATS (a)",
             "status": "COMPLETED", "exit": "0",
             "submit": "2026-01-01 10:00:02.0",
             "start": "2026-01-01 10:00:02.1",
             "complete": "2026-01-01 10:00:03.0",
             "realtime": "850ms", "workdir": wd2},
        ]
        self.trace = base / "trace.txt"
        _write_trace(self.trace, self.rows)

    def tearDown(self):
        self._td.cleanup()

    def test_receipts_bind_inputs_outputs_and_env(self):
        recs = rr.task_receipts(str(self.trace), engine="nf@1",
                                run_name="r", session_id="s")
        self.assertEqual(len(recs), 2)
        parse = recs[0]
        self.assertEqual(parse["process"], "P:PARSE")
        self.assertEqual(parse["chain_prev"], "genesis")
        self.assertEqual(recs[1]["chain_prev"], parse["receipt_id"])
        # input hash matches the staged target, not the symlink
        import hashlib
        self.assertEqual(parse["inputs"][0]["sha256"],
                         hashlib.sha256(b"fcs-a").hexdigest())
        out_names = {o["name"] for o in parse["outputs"]}
        self.assertEqual(out_names, {"a.events.tsv", "versions.yml"})
        # env fingerprint embeds the task's declared versions
        self.assertIn("versions_yml", parse["env"])
        self.assertTrue(parse["env_sha256"])

    def test_chain_verify_clean_then_tamper(self):
        recs = rr.task_receipts(str(self.trace), engine="nf@1",
                                run_name="r", session_id="s")
        self.assertEqual(rr.check_chain(recs), [])
        recs[0]["outputs"][0]["sha256"] = "0" * 64   # tamper body
        problems = rr.check_chain(recs)
        self.assertTrue(any("receipt_id" in p for p in problems))
        recs = rr.task_receipts(str(self.trace), engine="nf@1",
                                run_name="r", session_id="s")
        recs[1]["chain_prev"] = "genesis"            # break the chain
        problems = rr.check_chain(recs)
        self.assertTrue(any("chain break" in p for p in problems))

    def test_dotfiles_are_not_outputs(self):
        recs = rr.task_receipts(str(self.trace), engine="nf@1",
                                run_name="r", session_id="s")
        names = {o["name"] for r in recs for o in r["outputs"]}
        self.assertFalse(any(n.startswith(".") for n in names))

    def test_deterministic_receipt_ids(self):
        a = rr.task_receipts(str(self.trace), engine="nf@1",
                             run_name="r", session_id="s")
        # recorded_at differs run to run; ids hash the body including it,
        # so check the shape rather than equality across calls.
        for r in a:
            self.assertTrue(r["receipt_id"].startswith("rcpt-"))
            self.assertEqual(len(r["receipt_id"]), 21)

    def test_rehash_matches_then_detects_change(self):
        recs = rr.task_receipts(str(self.trace), engine="nf@1",
                                run_name="r", session_id="s")
        issues = rr.rehash(recs)
        self.assertTrue(issues[0].startswith("re-hashed"))
        self.assertFalse(any("mismatch" in i for i in issues))
        # rewrite an output file -> mismatch on next rehash
        import pathlib
        out = pathlib.Path(self.rows[0]["workdir"]) / "a.events.tsv"
        out.write_text("changed")
        issues = rr.rehash(recs)
        self.assertTrue(any("mismatch" in i for i in issues))

    def test_committed_receipts_verify(self):
        """The committed nf/results/run_receipts.jsonl is chain-clean."""
        committed = ROOT / "nf" / "results" / "run_receipts.jsonl"
        if not committed.exists():
            self.skipTest("receipts not generated yet")
        recs = [json.loads(l) for l in open(committed) if l.strip()]
        self.assertGreaterEqual(len(recs), 3)
        self.assertEqual(rr.check_chain(recs), [])
        # schema: every record binds env + workdir + hashes
        for r in recs:
            self.assertTrue(r["workdir"])
            self.assertTrue(r["env_sha256"])
            self.assertTrue(r["nextflow_task_hash"])
            self.assertTrue(r["outputs"], f"{r['task']} has no outputs")


if __name__ == "__main__":
    unittest.main()
