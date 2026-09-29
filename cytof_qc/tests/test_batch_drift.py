"""Batch alignment and acquisition-drift QC tests."""

import unittest

import numpy as np
import pandas as pd

from cytof_qc import batch, qc


def _frame(n=400, seed=0):
    rng = np.random.default_rng(seed)
    return pd.DataFrame(
        {"CD45": rng.normal(3.0, 0.2, n), "CD20": rng.normal(1.0, 0.2, n)}
    )


class AlignBatches(unittest.TestCase):
    def test_recovers_injected_offset(self):
        a = _frame(400, seed=1)
        b = _frame(400, seed=2)
        b["CD45"] += 0.8  # injected batch shift
        events = pd.concat(
            [a.assign(batch="H1"), b.assign(batch="H2")], ignore_index=True
        )
        out, report = batch.align_batches(events, ["CD45", "CD20"])
        # corrected batch medians should now sit near the pooled median
        for ch in ["CD45", "CD20"]:
            meds = out.groupby("batch")[ch].median()
            self.assertAlmostEqual(meds["H1"], meds["H2"], delta=0.05)
        # pooled median sits midway between batch medians, so the applied
        # shift is half the injected offset
        self.assertAlmostEqual(report["CD45"]["H2"]["shift"], -0.4, delta=0.05)
        self.assertAlmostEqual(report["CD45"]["H1"]["shift"], 0.4, delta=0.05)

    def test_tiny_batch_skipped_not_shifted(self):
        a = _frame(400)
        b = _frame(20, seed=3)
        events = pd.concat(
            [a.assign(batch="H1"), b.assign(batch="H2")], ignore_index=True
        )
        out, report = batch.align_batches(events, ["CD45"])
        self.assertTrue(report["CD45"]["H2"]["skipped"])
        # untouched: H2 values equal input
        np.testing.assert_allclose(
            out.loc[out.batch == "H2", "CD45"].to_numpy(),
            events.loc[events.batch == "H2", "CD45"].to_numpy(),
        )


class SplitChannels(unittest.TestCase):
    def test_qc_columns_excluded_from_phenotype(self):
        cols = ["CD45", "DNA1(Ir191)Di", "Viability(Pt195)Di", "Time",
                "event_number", "CD3"]
        pheno, qconly = qc.split_channels(cols)
        self.assertEqual(pheno, ["CD45", "CD3"])
        self.assertEqual(
            set(qconly),
            {"DNA1(Ir191)Di", "Viability(Pt195)Di", "Time", "event_number"},
        )


class AcquisitionQC(unittest.TestCase):
    def test_stable_channel_passes(self):
        df = _frame(1000)
        report = qc.acquisition_qc(df, ["CD45", "CD20"])
        self.assertFalse(report["CD45"]["flag_drift"])
        self.assertFalse(report["CD20"]["flag_drift"])

    def test_injected_drift_flags(self):
        df = _frame(1000)
        ramp = np.linspace(0, 2.0, len(df))  # decaying detector
        df["CD45"] = df["CD45"] + ramp
        report = qc.acquisition_qc(df, ["CD45", "CD20"])
        self.assertTrue(report["CD45"]["flag_drift"])
        self.assertFalse(report["CD20"]["flag_drift"])
        self.assertEqual(len(report["CD45"]["bin_medians"]), 10)

    def test_batch_seam_not_mistaken_for_drift(self):
        # two batches offset by a constant look like a single drift event
        # when pooled; per-batch measurement must keep them clean
        a = _frame(500, seed=1)
        b = _frame(500, seed=2)
        b["CD45"] += 1.5
        events = pd.concat(
            [a.assign(batch="H1"), b.assign(batch="H2")], ignore_index=True
        )
        pooled = qc.acquisition_qc(events, ["CD45"])
        per = qc.acquisition_qc(events, ["CD45"], batch_col="batch")
        self.assertTrue(pooled["CD45"]["flag_drift"])
        self.assertFalse(per["H1"]["CD45"]["flag_drift"])
        self.assertFalse(per["H2"]["CD45"]["flag_drift"])


if __name__ == "__main__":
    unittest.main()
