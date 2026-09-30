"""Shipped configs stay loadable and self-consistent."""

from pathlib import Path

import pytest
import yaml

CFG_DIR = Path(__file__).resolve().parent.parent / "config"


@pytest.mark.parametrize("name", ["config.yaml", "pbmc3k.yaml"])
def test_config_parses_with_required_keys(name):
    cfg = yaml.safe_load((CFG_DIR / name).read_text())
    for key in ("dataset", "qc", "embed", "paths"):
        assert key in cfg, f"{name} missing {key}"
    for key in ("min_genes_per_cell", "max_pct_mt", "min_cells_per_gene"):
        assert key in cfg["qc"]
    assert cfg["dataset"]["mode"] in ("demo", "pbmc3k", "h5ad")


def test_pbmc3k_writes_to_isolated_results_dir():
    demo = yaml.safe_load((CFG_DIR / "config.yaml").read_text())
    pbmc = yaml.safe_load((CFG_DIR / "pbmc3k.yaml").read_text())
    assert pbmc["dataset"]["mode"] == "pbmc3k"
    assert pbmc["paths"]["results_dir"] != demo["paths"]["results_dir"]
    assert pbmc["paths"]["data_dir"] != demo["paths"]["data_dir"]
