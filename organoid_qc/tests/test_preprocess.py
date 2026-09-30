import anndata as ad
import numpy as np
import pandas as pd
import pytest
import yaml

from organoid_qc.data.load import load
from organoid_qc.qc.preprocess import preprocess


def _cfg(**over):
    with open("config/config.yaml") as f:
        cfg = yaml.safe_load(f)
    for k, v in over.items():
        cfg[k].update(v) if isinstance(v, dict) else cfg.update({k: v})
    return cfg


def _mini_adata(n_obs=60, n_vars=80, cell_types=("a", "b"), seed=0):
    rng = np.random.default_rng(seed)
    X = rng.poisson(2.0, size=(n_obs, n_vars)).astype(np.float32)
    # a few high-signal genes so HVG has something to pick
    X[:, :10] *= 20
    obs = pd.DataFrame({
        "cell_type": np.tile(list(cell_types), n_obs // len(cell_types)
                             + 1)[:n_obs],
    }, index=[f"c{i}" for i in range(n_obs)])
    return ad.AnnData(X, obs=obs,
                      var=pd.DataFrame(index=[f"g{i}" for i in range(n_vars)]))


def _pp_cfg():
    # small but legal values for tiny fixtures
    return {"preprocess": {
        "min_genes_per_cell": 1, "min_cells_per_gene": 1,
        "target_sum": 10000, "n_hvg": 20, "n_pcs": 10,
        "leiden_resolution": 1.0}}


class TestLoadModes:
    def test_demo_mode_and_marker_carry(self):
        cfg = _cfg()
        organoid, reference = load(cfg)
        # demo mode keeps the synthetic panel (real-gene names in
        # cfg["markers"] would not exist in the synthetic matrix)
        assert organoid.uns["markers"] == cfg["dataset"]["demo"]["markers"]
        assert "cell_type" in reference.obs.columns

    def test_unknown_mode(self):
        cfg = _cfg(dataset={"mode": "bogus"})
        with pytest.raises(ValueError, match="dataset mode"):
            load(cfg)

    def test_missing_cell_type_column_rejected(self):
        cfg = _cfg(dataset={
            "mode": "h5ad",
            "organoid_h5ad": "x.h5ad",
            "reference_h5ad": "x.h5ad",
            "obs_columns": {"cell_type": "missing_col",
                            "cluster": "leiden"},
        })
        # h5ad mode reads files; point the read at a real in-memory
        # round trip by writing a small fixture instead
        import tempfile, os
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "r.h5ad")
            ref = _mini_adata()
            ref.obs = ref.obs.drop(columns=["cell_type"]) \
                          .assign(wrong="x")
            ref.write_h5ad(p)
            org = _mini_adata(seed=1)
            org.write_h5ad(os.path.join(d, "o.h5ad"))
            cfg["dataset"]["organoid_h5ad"] = os.path.join(d, "o.h5ad")
            cfg["dataset"]["reference_h5ad"] = p
            with pytest.raises(ValueError, match="missing_col"):
                load(cfg)


class TestPreprocess:
    def test_shared_gene_space_and_union_hvg(self):
        cfg = _pp_cfg()
        organoid = _mini_adata(seed=1)
        reference = _mini_adata(seed=2)
        # reference has 20 genes the organoid lacks and vice versa
        organoid = ad.concat(
            [organoid, _mini_adata(n_obs=60, n_vars=20, seed=3)],
            axis=1, join="outer", fill_value=0)
        organoid.var_names = ([f"g{i}" for i in range(80)]
                              + [f"org_only{i}" for i in range(20)])
        o, r = preprocess(organoid, reference, cfg)
        assert o.var_names.equals(r.var_names)
        assert o.n_vars == 80  # intersection, not union
        # union HVG flag: identical mask on both objects
        assert (o.var["highly_variable"] == r.var["highly_variable"]).all()
        assert o.var["highly_variable"].sum() > 0

    def test_lognorm_layer_and_clusters(self):
        cfg = _pp_cfg()
        o, r = preprocess(_mini_adata(seed=4), _mini_adata(seed=5), cfg)
        assert "lognorm" in o.layers and "lognorm" in r.layers
        assert "leiden" in o.obs.columns and "leiden" in r.obs.columns
        # scaled X is zero-centered but lognorm keeps raw normalized values
        assert o.layers["lognorm"].min() >= 0

    def test_disjoint_gene_names_raise(self):
        cfg = _pp_cfg()
        organoid = _mini_adata(seed=6)
        reference = _mini_adata(seed=7)
        reference.var_names = [f"ref{i}" for i in range(reference.n_vars)]
        with pytest.raises(ValueError, match="share no gene names"):
            preprocess(organoid, reference, cfg)
