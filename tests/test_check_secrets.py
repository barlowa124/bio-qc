"""The pre-commit secret scanner: pattern coverage and allowlist."""

from __future__ import annotations

import subprocess
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCANNER = ROOT / "scripts" / "check_secrets.py"

import sys
sys.path.insert(0, str(ROOT / "scripts"))
import check_secrets  # noqa: E402


class ScanTextTests(unittest.TestCase):
    def test_hf_token_detected(self):
        tok = "hf_" + "a" * 34
        hits = check_secrets.scan_text(f'"{tok}"', "nb.ipynb")
        self.assertEqual(len(hits), 1)
        self.assertIn("Hugging Face", hits[0][2])

    def test_redacted_marker_allowed(self):
        self.assertEqual(
            check_secrets.scan_text('"hf_REDACTED"', "nb.ipynb"), [])

    def test_placeholder_bodies_allowed_not_real_tokens(self):
        # all-x bodies are doc placeholders
        self.assertEqual(check_secrets.scan_text(
            '"hf_' + "x" * 34 + '"', "d.md"), [])
        # mixed-alnum bodies are credential-shaped: must flag
        tok = "hf_" + "aB3xZ9" * 5
        self.assertEqual(len(check_secrets.scan_text(
            f'"{tok}"', "d.md")), 1)

    def test_aws_example_key_allowed(self):
        # AWS's published documentation example is not a credential.
        self.assertEqual(check_secrets.scan_text(
            "AKIAIOSFODNN7EXAMPLE", "doc.md"), [])
        self.assertEqual(len(check_secrets.scan_text(
            "AKIA" + "Z9Y8X7W6V5U4T3S2", "k.py")), 1)

    def test_short_sk_strings_not_flagged(self):
        # CSS classes and short fragments are not keys.
        self.assertEqual(check_secrets.scan_text('class="sk-card"', "x.css"), [])

    def test_scanner_source_does_not_self_trigger(self):
        # the staged scanner + this test file must not flag themselves
        for path in (SCANNER, Path(__file__)):
            self.assertEqual(
                check_secrets.scan_text(
                    path.read_text(), str(path)), [],
                f"{path.name} would block its own commit")

    def test_new_token_families(self):
        self.assertTrue(check_secrets.scan_text(
            "pypi-" + "AgEIcHlwaS5vcmc" * 3, "c.toml"))
        self.assertTrue(check_secrets.scan_text(
            "glpat-" + "a1B2c3D4e5F6g7H8i9J0", "c.cfg"))
        self.assertTrue(check_secrets.scan_text(
            "npm_" + "aB" * 18, "c.ini"))

    def test_private_key_block_detected(self):
        # split literal so the scanner doesn't flag this test file
        hits = check_secrets.scan_text(
            "-----BEGIN " + "OPENSSH PRIVATE KEY-----", "id")
        self.assertTrue(hits)

    def test_clean_python_passes(self):
        self.assertEqual(check_secrets.scan_text(
            'token = os.environ["API_TOKEN"]\n', "m.py"), [])


class StagedModeTests(unittest.TestCase):
    def test_staged_scan_flags_added_token(self):
        import tempfile
        with tempfile.TemporaryDirectory() as td:
            subprocess.run(["git", "init", "-q"], cwd=td, check=True)
            subprocess.run(["git", "config", "user.email", "t@t"],
                           cwd=td, check=True)
            subprocess.run(["git", "config", "user.name", "t"],
                           cwd=td, check=True)
            Path(td, "clean.py").write_text("x = 1\n")
            subprocess.run(["git", "add", "."], cwd=td, check=True)
            subprocess.run(["git", "commit", "-qm", "init"],
                           cwd=td, check=True)
            Path(td, "nb.py").write_text(
                f'key = "hf_{"z" * 34}"\nok = "fine"\n')
            subprocess.run(["git", "add", "nb.py"], cwd=td, check=True)
            r = subprocess.run(
                [sys.executable, str(SCANNER), "--staged"],
                cwd=td, capture_output=True, text=True)
            self.assertEqual(r.returncode, 1)
            self.assertIn("nb.py", r.stderr)

    def test_staged_scan_ignores_clean_tree(self):
        import tempfile
        with tempfile.TemporaryDirectory() as td:
            subprocess.run(["git", "init", "-q"], cwd=td, check=True)
            r = subprocess.run(
                [sys.executable, str(SCANNER), "--staged"],
                cwd=td, capture_output=True, text=True)
            self.assertEqual(r.returncode, 0)


if __name__ == "__main__":
    unittest.main()
