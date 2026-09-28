"""The vendored canonical conformal module vs the notebook's inline
quantile formula.

notebooks/conformal_prediction.py uses np.quantile with default linear
interpolation at the ceiled level, and raises when the level exceeds 1
(small n, small alpha). The shared module reads the order statistic
directly and stays defined for all n.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import conformal


def test_canonical_never_below_notebook_formula():
    rng = np.random.default_rng(0)
    for _ in range(200):
        n = int(rng.integers(20, 400))
        s = rng.exponential(1.0, n)
        alpha = float(rng.choice([0.05, 0.1, 0.2]))
        level = np.ceil((1 + n) * (1 - alpha)) / n
        notebook_q = min(np.quantile(s, level), 1.0)
        assert conformal.qhat(s, alpha) >= notebook_q - 1e-12


def test_small_n_small_alpha_defined():
    # The notebook formula raises at level>1 (n<19, alpha=0.05); the
    # canonical level is the largest order statistic.
    s = np.array([0.1, 0.2, 0.3, 0.4, 0.5])
    assert conformal.qhat(s, 0.05) == 0.5


def test_class_helpers():
    probs = np.array([[0.7, 0.3], [0.4, 0.6]])
    y = np.array([0, 1])
    np.testing.assert_allclose(conformal.class_scores(probs, y),
                               [0.3, 0.4])
    assert list(conformal.covered_set(y, conformal.class_sets(probs, 0.4))) == [True, True]
