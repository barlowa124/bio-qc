"""Deterministic demo generator: structure, reproducibility, injected
low-quality cells actually look low-quality."""

import numpy as np
import pytest

from scrna_qc.data import make_demo


@pytest.fixture()
def cfg():
    return {"dataset": {"demo": {
        "seed": 7, "n_cells": 300, "n_genes": 400, "n_mt_genes": 12,
        "type_proportions": [0.5, 0.3, 0.2],
        "marker_block_size": 20,
        "frac_low_quality": 0.10, "frac_mito_high": 0.08}}}


def test_demo_deterministic(cfg):
    a = make_demo(cfg)
    b = make_demo(cfg)
    assert (a.X != b.X).nnz == 0
    assert a.obs["true_type"].tolist() == b.obs["true_type"].tolist()


def test_demo_structure(cfg):
    a = make_demo(cfg)
    assert a.n_obs == 300
    assert a.n_vars == 400
    assert sum(v.startswith("MT-") for v in a.var_names) == 12
    assert a.obs["true_type"].nunique() == 3


def test_demo_marker_blocks_distinct(cfg):
    a = make_demo(cfg)
    types = np.asarray(a.obs["true_type"])
    x0 = np.asarray(a.X[types == "type_0"].todense()).mean(axis=0)
    x1 = np.asarray(a.X[types == "type_1"].todense()).mean(axis=0)
    # type_0's marker block sits at var index 12..32; type_1's at 32..52
    assert x0[12:32].mean() > 3 * x1[12:32].mean()


def test_injected_low_quality_cells_have_low_depth(cfg):
    a = make_demo(cfg)
    totals = np.asarray(a.X.sum(axis=1)).ravel()
    deg = np.asarray(a.obs["true_degraded"])
    healthy_median = np.median(totals[deg == "none"])
    assert np.all(totals[deg == "low_genes"] < 0.3 * healthy_median)


def test_injected_mito_cells_have_high_mt_fraction(cfg):
    from scrna_qc.qc_metrics import add_qc_metrics
    a = add_qc_metrics(make_demo(cfg), "MT-")
    deg = np.asarray(a.obs["true_degraded"])
    assert np.all(np.asarray(a.obs["pct_mt"])[deg == "high_mito"] > 25)
