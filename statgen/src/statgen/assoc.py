"""Association stage: single-variant tests with covariates.

linear   - OLS per variant via Frisch-Waugh-Lovell: residualize the
           phenotype and the whole genotype matrix on [1, cov, PCs]
           once, then beta/se/t for every variant from vectorized dot
           products. p from the t distribution.
logistic - Rao score test per variant against a null logistic fit
           (IRLS on [1, cov, PCs] only), chi-square 1 df.

Genomic-control lambda = median(chi2) / 0.4549364 is reported both
with covariates (the headline) and with no covariates, so the effect
of stratification control is measurable rather than asserted.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2, norm, t as tdist

from statgen.util import load_config

CHI2_MEDIAN = 0.4549364


def _impute(G: np.ndarray) -> np.ndarray:
    """Mean-impute missing calls per variant."""
    Gf = np.where(G < 0, np.nan, G).astype(np.float64)
    means = np.nanmean(Gf, axis=0)
    idx = np.where(np.isnan(Gf))
    Gf[idx] = np.take(means, idx[1])
    return Gf


def _residualize(M: np.ndarray, X: np.ndarray) -> np.ndarray:
    """Residuals of each column of M on X (least squares)."""
    beta, *_ = np.linalg.lstsq(X, M, rcond=None)
    return M - X @ beta


def linear_scan(y: np.ndarray, Gf: np.ndarray, X: np.ndarray):
    """Returns (beta, se, chi2_stat, p) per variant."""
    ry = _residualize(y, X)
    rG = _residualize(Gf, X)
    n, k = X.shape
    df = n - k - 1
    num = rG * ry[:, None]
    denom = np.sum(rG ** 2, axis=0)
    denom[denom == 0] = np.nan
    beta = num.sum(axis=0) / denom
    resid = ry[:, None] - rG * beta[None, :]
    sigma2 = np.sum(resid ** 2, axis=0) / df
    se = np.sqrt(sigma2 / denom)
    with np.errstate(divide="ignore", invalid="ignore"):
        tstat = beta / se
    stat = tstat ** 2
    p = 2.0 * tdist.sf(np.abs(tstat), df)
    return beta, se, stat, p


def _fit_null_logistic(y: np.ndarray, X: np.ndarray, iters: int = 25):
    b = np.zeros(X.shape[1])
    for _ in range(iters):
        eta = np.clip(X @ b, -30, 30)
        mu = 1.0 / (1.0 + np.exp(-eta))
        w = np.maximum(mu * (1 - mu), 1e-9)
        grad = X.T @ (y - mu)
        h = (X * w[:, None]).T @ X
        try:
            step = np.linalg.solve(h, grad)
        except np.linalg.LinAlgError:
            step = np.linalg.lstsq(h, grad, rcond=None)[0]
        b += step
        if np.max(np.abs(step)) < 1e-8:
            break
    eta = np.clip(X @ b, -30, 30)
    mu = 1.0 / (1.0 + np.exp(-eta))
    w = np.maximum(mu * (1 - mu), 1e-9)
    return mu, w


def logistic_score_scan(y: np.ndarray, Gf: np.ndarray, X: np.ndarray):
    """Rao score test for adding each variant to the null model."""
    mu, w = _fit_null_logistic(y, X)
    r = y - mu
    Xw = X * w[:, None]
    # efficient score: G'r - G'W X (X'WX)^-1 X'r
    XtWX_inv = np.linalg.pinv(Xw.T @ X)
    score_g = Gf.T @ r
    cross = Gf.T @ Xw
    info = (Gf * w[:, None]).T @ Gf - cross @ XtWX_inv @ cross.T
    var = np.diag(info).copy()
    var[var <= 0] = np.nan
    z = score_g / np.sqrt(var)
    stat = z ** 2
    # sign of the score gives the effect direction; |beta| approx =
    # score / info, standard logistic score-test reporting
    beta = score_g / var
    se = 1.0 / np.sqrt(var)
    p = chi2.sf(stat, 1)
    return beta, se, stat, p


def genomic_lambda(stat: np.ndarray) -> float:
    s = stat[np.isfinite(stat)]
    if s.size == 0:
        return float("nan")
    return float(np.median(s) / CHI2_MEDIAN)


def main() -> None:
    cfg = load_config()
    data_dir = Path(cfg["paths"]["data_dir"])
    results_dir = Path(cfg["paths"]["results_dir"])
    z = np.load(data_dir / "processed" / "filtered.npz",
                allow_pickle=True)
    pcs = pd.read_csv(results_dir / "pcs.tsv", sep="\t")
    a = cfg["assoc"]
    n_pcs = int(a["n_pcs_as_cov"])

    y = z["y"].astype(float)
    Gf = _impute(z["G"])
    cov = z["cov"].astype(float) if "cov" in z else np.zeros((len(y), 0))
    X = np.column_stack(
        [np.ones(len(y)), cov,
         pcs[[f"PC{i+1}" for i in range(n_pcs)]].to_numpy()])
    X0 = np.ones((len(y), 1))

    kind = a["kind"]
    y_eval = y
    if kind == "linear":
        scan = linear_scan
    elif kind == "logistic":
        y_eval = (y > np.median(y)).astype(float)
        scan = logistic_score_scan
    else:
        raise ValueError(f"unknown assoc.kind: {kind}")

    beta, se, stat, p = scan(y_eval, Gf, X)
    _, _, stat0, _ = scan(y_eval, Gf, X0)

    # Dose-response curve: lambda at 0..n_pcs covariate PCs shows
    # how much stratification each PC absorbs, not just endpoints.
    pc_cols = pcs[[f"PC{i+1}" for i in range(n_pcs)]].to_numpy()
    lam_curve = {}
    for k in range(n_pcs + 1):
        Xk = np.column_stack([np.ones(len(y)), cov, pc_cols[:, :k]])
        _, _, stat_k, _ = scan(y_eval, Gf, Xk)
        lam_curve[f"pcs_{k}"] = genomic_lambda(stat_k)

    df = pd.DataFrame({
        "vid": z["vid"].astype(str), "chrom": z["chrom"].astype(str),
        "pos": z["pos"].astype(int),
        "a1": z["a1"].astype(str), "a2": z["a2"].astype(str),
        "beta": beta, "se": se, "chi2": stat, "p": p,
    })
    results_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(results_dir / "assoc.tsv", sep="\t", index=False)

    bonf = 0.05 / len(df)
    lead = df.loc[df["p"].idxmin()]
    summary = {
        "kind": kind,
        "n_variants": int(len(df)),
        "n_samples": int(len(y)),
        "n_covariates": int(X.shape[1] - 1),
        "n_pcs_as_cov": n_pcs,
        "lambda_gc": genomic_lambda(stat),
        "lambda_gc_no_cov": genomic_lambda(stat0),
        "lambda_by_pcs": lam_curve,
        "bonferroni_alpha": float(bonf),
        "n_bonferroni_hits": int((p < bonf).sum()),
        "lead_hit": {
            "vid": str(lead["vid"]), "pos": int(lead["pos"]),
            "chrom": int(lead["chrom"])
                if str(lead["chrom"]).isdigit() else str(lead["chrom"]),
            "beta": float(lead["beta"]), "se": float(lead["se"]),
            "p": float(lead["p"]),
        },
    }
    with open(results_dir / "assoc_summary.json", "w") as f:
        json.dump(summary, f, indent=2)


if __name__ == "__main__":
    sys.exit(main())
