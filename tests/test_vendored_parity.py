"""Vendored parity: conformal.py copies share the portfolio digest.

The canonical conformal helpers are vendored (not depended on) in
protein-ml/protein_stability_uncertainty, mol-ml/comp_tox_pipeline,
cultivated_meat_multiomic, and scrna_qc — four consumers. mol-ml and
protein-ml pin the same sha256; an edit anywhere trips a pin and forces
deliberate sync.
"""
import hashlib
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

COPIES = [
    ROOT / "cultivated_meat_multiomic" / "conformal.py",
    ROOT / "scrna_qc" / "src" / "scrna_qc" / "conformal.py",
]

# Shared with mol-ml and protein-ml tests — keep in sync.
CONFORMAL_SHA256 = \
    "fcba54721a3f106e3864ab43ad0cb61caf273c41426d7bfa649b542f5f72af00"


CLAIMS_COPIES = [
    ROOT / "scrna_qc" / "src" / "scrna_qc" / "claims.py",
    ROOT / "statgen" / "src" / "statgen" / "claims.py",
]

# sha256 of the claims verifier as ported from oncology_coscientist —
# bio-qc-internal pin; two consumers, drift fails here.
CLAIMS_SHA256 = \
    "475eb4a6a338e363374c0810bf8f3861aec16cd41ae0c4ec2e0fbb54a6cc74a2"


class VendoredParityTests(unittest.TestCase):
    def test_conformal_matches_shared_digest(self):
        for path in COPIES:
            self.assertEqual(
                hashlib.sha256(path.read_bytes()).hexdigest(),
                CONFORMAL_SHA256,
                f"{path.relative_to(ROOT)} drifted from the vendored "
                "copies — sync all copies and update the pin together")

    def test_claims_matches_shared_digest(self):
        for path in CLAIMS_COPIES:
            self.assertEqual(
                hashlib.sha256(path.read_bytes()).hexdigest(),
                CLAIMS_SHA256,
                f"{path.relative_to(ROOT)} drifted from the vendored "
                "claims verifier — sync scrna_qc and statgen together")


if __name__ == "__main__":
    unittest.main()
