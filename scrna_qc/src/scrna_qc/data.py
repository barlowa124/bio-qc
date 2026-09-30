"""Dataset loading.

- demo: deterministic synthetic AnnData. Four populations with uneven
  proportions and distinct marker blocks, plus injected low-quality
  (few detected genes) and high-mito cells so the filters have real work
  to do. Synthetic by design; never a biological claim.
- pbmc3k: the public 10x PBMC 3k dataset via scanpy's cache.
- h5ad: a local file path from config.
"""

from __future__ import annotations

import sys
from pathlib import Path

import anndata as ad
import numpy as np
from scipy import sparse

from .util import load_config

MT_FRACTION_STRESSED = 0.40
LOW_QUALITY_SCALE = 0.02


def make_demo(cfg: dict) -> ad.AnnData:
    d = cfg["dataset"]["demo"]
    rng = np.random.default_rng(d["seed"])
    n_genes = d["n_genes"]
    n_mt = d["n_mt_genes"]

    var_names = [f"MT-{i + 1:04d}" for i in range(n_mt)] + [
        f"G{i + 1:04d}" for i in range(n_genes - n_mt)
    ]

    base = rng.gamma(shape=1.2, scale=0.6, size=n_genes)
    base[:n_mt] *= 1.5  # mt genes at a moderate baseline

    type_props = d["type_proportions"]
    block = d["marker_block_size"]
    type_profiles = []
    for t, _ in enumerate(type_props):
        p = base.copy()
        start = n_mt + t * block
        p[start:start + block] *= 8.0
        type_profiles.append(p / p.sum())

    rows = []
    true_type = []
    degraded = []
    for t, frac in enumerate(type_props):
        n_t = round(d["n_cells"] * frac)
        size = rng.lognormal(mean=8.0, sigma=0.35, size=n_t)
        lam = size[:, None] * type_profiles[t][None, :]
        x = rng.poisson(lam)
        true_type += [f"type_{t}"] * n_t
        degraded += ["none"] * n_t
        rows.append(x)
    X = np.vstack(rows).astype(np.float64)

    n_cells = X.shape[0]
    n_low = round(n_cells * d["frac_low_quality"])
    n_mito = round(n_cells * d["frac_mito_high"])
    idx = rng.permutation(n_cells)

    for i in idx[:n_low]:
        X[i] = rng.poisson(X[i] * LOW_QUALITY_SCALE)
        degraded[i] = "low_genes"
    for i in idx[n_low:n_low + n_mito]:
        non_mt = X[i, n_mt:]
        mt_total = X[i].sum() * MT_FRACTION_STRESSED / (1 - MT_FRACTION_STRESSED)
        X[i, :n_mt] = rng.poisson(
            np.full(n_mt, mt_total / n_mt)
        )
        X[i, n_mt:] = non_mt * 0.6
        degraded[i] = "high_mito"

    adata = ad.AnnData(
        X=sparse.csr_matrix(X),
        obs={"true_type": true_type, "true_degraded": degraded},
        var={"gene": var_names},
    )
    adata.obs_names = [f"CELL{i + 1:05d}" for i in range(n_cells)]
    adata.var_names = var_names
    return adata


def load_dataset(cfg: dict) -> ad.AnnData:
    mode = cfg["dataset"]["mode"]
    if mode == "demo":
        return make_demo(cfg)
    if mode == "pbmc3k":
        import scanpy as sc
        return sc.datasets.pbmc3k()
    if mode == "h5ad":
        path = cfg["dataset"]["h5ad_path"]
        if not path:
            raise ValueError("dataset.mode is h5ad but h5ad_path is null")
        return ad.read_h5ad(path)
    raise ValueError(f"unknown dataset mode: {mode}")


def main() -> None:
    out = Path(sys.argv[1])
    out.parent.mkdir(parents=True, exist_ok=True)
    adata = load_dataset(load_config())
    adata.write_h5ad(out)
    print(f"data: {adata.n_obs} cells x {adata.n_vars} genes -> {out}")


if __name__ == "__main__":
    main()
