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


def pca(Z: np.ndarray, k: int) -> tuple[np.ndarray, np.ndarray,
                                       np.ndarray]:
    """Top-k sample coords, eigenvalues of Z Z^T / (V-1), and the
    variant loadings (right singular vectors)."""
    v = Z.shape[1]
    _, s, vt = np.linalg.svd(Z, full_matrices=False)
    eig = (s ** 2) / max(v - 1, 1)
    pcs = Z @ vt[:k].T
    return pcs, eig, vt


def kinship_pairs(Z: np.ndarray, loadings: np.ndarray, k_proj: int,
                  flag: float) -> tuple[float, int]:
    """Genomic relatedness off-diagonal scan on PC-projected
    genotypes (the PC-AiR idea, simplified): in a structured cohort the
    raw Z Z^T / V off-diagonal captures ancestry sharing, not recent
    relatedness, so the top variant-loadings are projected out first.
    Pairs above `flag` (0.125 ~ third-degree) are cryptic relatives.
    Returns (max off-diagonal, off-diagonal sd, flagged pair count).
    The sd matters: with few or LD-dense variants the estimate is
    overdispersed and a high flag count means the estimate is
    unreliable, not that the cohort is full of relatives."""
    n, v = Z.shape
    vk = loadings[:k_proj]
    Zr = Z - (Z @ vk.T) @ vk
    grm = (Zr @ Zr.T) / v
    iu = np.triu_indices(n, k=1)
    off = grm[iu]
    return float(off.max()), float(off.std()), int((off > flag).sum())


def run(G: np.ndarray, n_pcs: int) -> tuple[np.ndarray, np.ndarray]:
    pcs, eig, _ = pca(standardize(G), n_pcs)
    return pcs, eig


def main() -> None:
    cfg = load_config()
    data_dir = Path(cfg["paths"]["data_dir"])
    results_dir = Path(cfg["paths"]["results_dir"])
    z = np.load(data_dir / "processed" / "filtered.npz",
                allow_pickle=True)
    n_pcs = int(cfg["strat"]["n_pcs"])

    Z = standardize(z["G"])
    pcs, eig, vt = pca(Z, n_pcs)
    var_expl = eig / eig.sum()
    k_flag = float(cfg["strat"].get("kinship_flag", 0.125))
    kin_max, kin_sd, n_related = kinship_pairs(Z, vt, n_pcs, k_flag)

    summary = {
        "n_pcs": n_pcs,
        "var_explained": [float(x) for x in var_expl[:n_pcs]],
        "top_eigenvalues": [float(x) for x in eig[:n_pcs]],
        "n_variants_used": int(z["G"].shape[1]),
        "n_samples": int(z["G"].shape[0]),
        "kinship_flag": k_flag,
        "kinship_offdiag_max": kin_max,
        "kinship_offdiag_sd": kin_sd,
        "n_related_pairs": n_related,
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
