"""split_audit: planted-leakage recovery and clean-split acceptance."""
import unittest

from split_audit import audit, demo


class ScaffoldAuditTests(unittest.TestCase):
    def test_demo_planted_scaffolds_recovered(self):
        scaf_rows, _, _ = demo.make_demo()
        result = audit.audit_scaffolds(scaf_rows)
        self.assertEqual(2, result["n_shared_scaffolds"])
        self.assertFalse(result["scaffold_disjoint"])
        # each shared scaffold contributes its train row plus its test row
        self.assertEqual(4, result["n_leaked_items"])

    def test_clean_split_is_disjoint(self):
        rows = [{"id": f"m{i}", "scaffold_id": f"s{i}",
                 "split": "train" if i < 5 else "test"}
                for i in range(10)]
        result = audit.audit_scaffolds(rows)
        self.assertTrue(result["scaffold_disjoint"])
        self.assertEqual(0, result["n_shared_scaffolds"])
        self.assertEqual(0, result["n_leaked_items"])


class SequenceAuditTests(unittest.TestCase):
    def test_demo_planted_sequences_recovered(self):
        _, seq_rows, planted = demo.make_demo()
        result = audit.audit_sequences(seq_rows)
        flagged = {p["a"] for p in result["flagged_pairs"]} \
            | {p["b"] for p in result["flagged_pairs"]}
        self.assertTrue(set(planted) <= flagged)
        self.assertFalse(result["leakage_free"])

    def test_clean_sequences_pass(self):
        rng_rows = [{"id": "a", "sequence": "A" * 30, "split": "train"},
                    {"id": "b", "sequence": "C" * 30, "split": "test"}]
        result = audit.audit_sequences(rng_rows)
        self.assertTrue(result["leakage_free"])
        self.assertEqual(0, result["n_flagged_pairs"])

    def test_threshold_boundary(self):
        rows = [{"id": "x", "sequence": "ACDE" * 20, "split": "train"},
                {"id": "y", "sequence": "ACDE" * 20, "split": "test"}]
        self.assertEqual(1, audit.audit_sequences(
            rows, threshold=0.999)["n_flagged_pairs"])
        self.assertEqual(0, audit.audit_sequences(
            rows, threshold=1.0001)["n_flagged_pairs"])


if __name__ == "__main__":
    unittest.main()
