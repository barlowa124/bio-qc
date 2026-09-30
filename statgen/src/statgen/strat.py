"""Population stratification: PCA on the standardized genotype matrix.

Standardization follows the EIGENSOFT convention: each variant is
centered on its mean dosage and scaled by sqrt(2 p (1-p)); missing
calls are imputed to the variant mean (zero after centering).

Outputs:
  pcs.tsv         S x K sample coordinates (used downstream as
                  association covariates)
  strat_summary.json  eigenvalues, variance explained, and — on the
                  demo cohort where ancestry labels exist — the
                  |r| between PC1 and the ancestry label so the
                  structure being corrected is measurable.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from statgen.util import load_config


def standardize(G: np.ndarray) -> np.ndarray:
    Gf = np.where(G < 0, np.nan, G).astype(np.float64)
    p = np.nanmean(Gf, axis=0) / 2.0
    scale = np.sqrt(2.0 * p * (1.0 - p))
    scale[scale == 0] = 1.0
    Z = (np.nan_to_num(Gf, nan=np.nan) - 2.0 * p) / scale
    Z[np.isnan(Z)] = 0.0
    return Z


def pca(Z: np.ndarray, k: int) -> tuple[np.ndarray, np.ndarray]:
    """Top-k right-singular coords and eigenvalues of Z Z^T / (V-1)."""
    v = Z.shape[1]
    _, s, vt = np.linalg.svd(Z, full_matrices=False)
    eig = (s ** 2) / max(v - 1, 1)
    pcs = Z @ vt[:k].T
    return pcs, eig


def run(G: np.ndarray, n_pcs: int) -> tuple[np.ndarray, np.ndarray]:
    return pca(standardize(G), n_pcs)


def main() -> None:
    cfg = load_config()
    data_dir = Path(cfg["paths"]["data_dir"])
    results_dir = Path(cfg["paths"]["results_dir"])
    z = np.load(data_dir / "processed" / "filtered.npz",
                allow_pickle=True)
    n_pcs = int(cfg["strat"]["n_pcs"])

    pcs, eig = run(z["G"], n_pcs)
    var_expl = eig / eig.sum()

    summary = {
        "n_pcs": n_pcs,
        "var_explained": [float(x) for x in var_expl[:n_pcs]],
        "top_eigenvalues": [float(x) for x in eig[:n_pcs]],
        "n_variants_used": int(z["G"].shape[1]),
        "n_samples": int(z["G"].shape[0]),
    }
    if "ancestry" in z and z["ancestry"].size:
        anc = z["ancestry"].astype(float)
        groups = np.unique(anc)
        if len(groups) == 2:
            r = np.corrcoef(pcs[:, 0], anc)[0, 1]
            summary["pc1_ancestry_r"] = float(r)
        # eta^2: share of PC1 variance between groups (works for
        # arbitrary label orderings and multi-group panels)
        tot = np.var(pcs[:, 0])
        means = np.array([pcs[anc == a, 0].mean() for a in groups])
        ns = np.array([(anc == a).sum() for a in groups])
        between = float((ns * (means - pcs[:, 0].mean()) ** 2).sum()
                        / len(anc))
        summary["pc1_ancestry_eta2"] = (
            between / tot if tot > 0 else 0.0)
        summary["n_ancestry_groups"] = int(len(groups))

    results_dir.mkdir(parents=True, exist_ok=True)
    cols = {f"PC{i+1}": pcs[:, i] for i in range(n_pcs)}
    df = pd.DataFrame({"sample": z["samples"].astype(str), **cols})
    df.to_csv(results_dir / "pcs.tsv", sep="\t", index=False)
    with open(results_dir / "strat_summary.json", "w") as f:
        json.dump(summary, f, indent=2)


if __name__ == "__main__":
    sys.exit(main())
