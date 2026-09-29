"""Load or generate the organoid + reference pair.

Modes (config.dataset.mode):
  demo  - deterministic synthetic data with known ground truth
  h5ad  - user-supplied annotated AnnData files (e.g. CELLxGENE Census
          exports); paths in config.dataset.organoid_h5ad / reference_h5ad

Writes both to data/processed/*.h5ad so downstream stages are mode-agnostic.
"""

from __future__ import annotations

import sys

import anndata as ad
from organoid_qc.util import load_config


def load(cfg: dict) -> tuple[ad.AnnData, ad.AnnData]:
    mode = cfg["dataset"]["mode"]
    if mode == "demo":
        from organoid_qc.data.synthetic import make_demo

        organoid, reference = make_demo(cfg)
    elif mode == "fetch":
        from organoid_qc.data.fetch import fetch

        organoid, reference = fetch(cfg)
    elif mode == "h5ad":
        organoid = ad.read_h5ad(cfg["dataset"]["organoid_h5ad"])
        reference = ad.read_h5ad(cfg["dataset"]["reference_h5ad"])
    else:
        raise ValueError(f"unknown dataset mode: {mode!r}")

    # Carry the active marker panel so scoring is mode-agnostic.
    organoid.uns.setdefault("markers", dict(cfg["markers"]))
    reference.uns.setdefault("markers", dict(cfg["markers"]))

    ct_col = cfg["dataset"]["obs_columns"]["cell_type"]
    if ct_col not in reference.obs.columns:
        raise ValueError(
            f"reference is missing required obs column {ct_col!r}"
        )
    return organoid, reference


def main() -> None:
    organoid_out, reference_out = sys.argv[1], sys.argv[2]
    cfg = load_config()
    organoid, reference = load(cfg)
    organoid.write_h5ad(organoid_out)
    reference.write_h5ad(reference_out)
    print(
        f"load: organoid {organoid.n_obs}x{organoid.n_vars}, "
        f"reference {reference.n_obs}x{reference.n_vars}"
    )


if __name__ == "__main__":
    main()
