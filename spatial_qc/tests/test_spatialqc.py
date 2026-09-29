import os
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import numpy as np  # noqa: E402

import spatialqc  # noqa: E402


def _write_sample(d: str) -> None:
    """Minimal Space Ranger export: 3 genes x 4 spots."""
    mm = Path(d) / "filtered_feature_bc_matrix"
    mm.mkdir(parents=True)
    (mm / "barcodes.tsv").write_text("BC1\nBC2\nBC3\nBC4\n")
    (mm / "features.tsv").write_text(
        "ENSG1\tGeneA\tGene Expression\n"
        "ENSG2\tmt-Nd1\tGene Expression\n"
        "ENSG3\tGeneC\tGene Expression\n")
    # genes x barcodes: spot1 has all genes, spot4 has 1 count 1 gene
    (mm / "matrix.mtx").write_text(
        "%%MatrixMarket matrix coordinate integer general\n"
        "3 4 7\n"
        "1 1 100\n1 2 80\n2 1 10\n2 3 60\n"
        "3 1 90\n3 2 70\n3 4 1\n")
    (Path(d) / "spatial").mkdir(exist_ok=True)
    (Path(d) / "spatial" / "tissue_positions_list.csv").write_text(
        "BC1,1,0,0,100,100\nBC2,1,0,1,100,200\n"
        "BC3,0,1,0,300,100\nBC4,1,1,1,300,200\n")


class LoaderTests(unittest.TestCase):
    def test_load(self):
        with tempfile.TemporaryDirectory() as d:
            _write_sample(d)
            s = spatialqc.load_visium(d)
        self.assertEqual(s["matrix"].shape, (4, 3))
        self.assertEqual(s["genes"], ["GeneA", "mt-Nd1", "GeneC"])
        self.assertEqual(s["positions"]["BC1"]["in_tissue"], "1")


class MetricsTests(unittest.TestCase):
    def test_metrics(self):
        with tempfile.TemporaryDirectory() as d:
            _write_sample(d)
            s = spatialqc.load_visium(d)
        met = spatialqc.spot_metrics(s)
        # BC1: 100+10+90=200 counts, 3 genes, mito frac 10/200
        self.assertEqual(met["total_counts"][0], 200)
        self.assertEqual(met["n_genes"][0], 3)
        self.assertAlmostEqual(met["mito_frac"][0], 0.05)
        # BC4: 1 count, 1 gene
        self.assertEqual(met["total_counts"][3], 1)
        self.assertEqual(met["in_tissue"].tolist(), [1, 1, 0, 1])
        summ = spatialqc.summary(s, met)
        self.assertEqual(summ["n_spots"], 4)
        self.assertEqual(summ["n_in_tissue"], 3)


class FilterTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        _write_sample(self._tmp.name)
        self.sample = spatialqc.load_visium(self._tmp.name)
        self.met = spatialqc.spot_metrics(self.sample)

    def tearDown(self):
        self._tmp.cleanup()

    def test_fixed_drops_low_spot(self):
        # default cutoffs (500 counts / 200 genes) drop all 4 tiny spots
        self.assertFalse(spatialqc.fixed(self.met).any())
        # spot-appropriate cutoffs still drop BC3 (60 counts but 100%
        # mitochondrial) and BC4 (1 count, 1 gene)
        keep = spatialqc.fixed(self.met, min_counts=2, min_genes=1)
        self.assertEqual(keep.tolist(), [True, True, False, False])

    def test_mad_adapts(self):
        keep = spatialqc.mad(self.met, nmad=5)
        self.assertEqual(len(keep), 4)
        # BC4 (1 count, 1 gene) is the distribution outlier
        self.assertFalse(keep[3])

    def test_tissue_only(self):
        keep = spatialqc.tissue_only(self.met)
        self.assertEqual(keep.tolist(), [True, True, False, True])

    def test_compare_table(self):
        rep = spatialqc.compare(self.sample, run_clusters=False)
        self.assertEqual(rep["n_spots_total"], 4)
        for name in ("fixed", "mad", "tissue_only"):
            self.assertIn(name, rep["strategies"])
        self.assertIn("spots_kept", rep["strategies"]["fixed"])
        self.assertIn("counts_kept_frac",
                      rep["strategies"]["mad"])


class CliTests(unittest.TestCase):
    def test_cli(self):
        from spatialqc.cli import main
        with tempfile.TemporaryDirectory() as d:
            _write_sample(d)
            out = os.path.join(d, "r.json")
            rc = main([d, "--out", out, "--no-cluster"])
            self.assertEqual(rc, 0)
            import json
            rep = json.load(open(out))
            self.assertEqual(rep["summary"]["n_spots"], 4)


if __name__ == "__main__":
    unittest.main()
