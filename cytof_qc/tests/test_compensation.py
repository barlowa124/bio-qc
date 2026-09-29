"""Compensation module tests."""

import unittest

import numpy as np
import pandas as pd

from cytof_qc import compensation


class Compensate(unittest.TestCase):
    def test_identity_is_noop(self):
        df = pd.DataFrame({"CD45": [1.0, 2.0], "CD20": [0.5, -0.1]})
        out = compensation.compensate(
            df, ["CD45", "CD20"], np.eye(2)
        )
        np.testing.assert_allclose(
            out[["CD45", "CD20"]].to_numpy(), df.to_numpy(), atol=1e-12
        )

    def test_recovers_known_spill(self):
        # 10% of CD45 observed in CD19; true signal must be recovered
        S = np.array([[1.0, 0.1], [0.0, 1.0]])
        true = np.array([[3.0, 0.0], [5.0, 0.2]])
        observed = true @ S
        df = pd.DataFrame(observed, columns=["CD45", "CD19"])
        out = compensation.compensate(df, ["CD45", "CD19"], S)
        np.testing.assert_allclose(out.to_numpy(), true, atol=1e-10)

    def test_wrong_shape_rejected(self):
        df = pd.DataFrame({"A": [1.0], "B": [2.0]})
        with self.assertRaises(ValueError):
            compensation.compensate(df, ["A", "B"], np.eye(3))

    def test_singular_rejected(self):
        df = pd.DataFrame({"A": [1.0], "B": [2.0]})
        with self.assertRaises(ValueError):
            compensation.compensate(
                df, ["A", "B"], np.array([[1.0, 1.0], [1.0, 1.0]])
            )


if __name__ == "__main__":
    unittest.main()
