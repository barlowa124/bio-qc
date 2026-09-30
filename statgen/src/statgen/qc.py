"""Genotype QC stage: per-sample and per-variant filters, all thresholds
from config. Writes the filtered cohort plus a reconciling waterfall:

  variants_in = sum over steps of variants_removed + variants_out
  samples_in  = sum over steps of samples_removed  + samples_out
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import chi2

from statgen.util import load_config


def maf(G: np.ndarray) -> np.ndarray:
    """Minor allele frequency per variant over called genotypes."""
    Gf = G.astype(float)
    Gf[Gf < 0] = np.nan
    with np.errstate(invalid="ignore"):
        freq = np.nanmean(Gf, axis=0) / 2.0
    return np.minimum(freq, 1.0 - freq)


def hwe_p(G: np.ndarray) -> np.ndarray:
    """Chi-square HWE p-value per variant (1 df), NaN on monomorphic."""
    n_v = G.shape[1]
    out = np.full(n_v, np.nan)
    for j in range(n_v):
        g = G[:, j]
        g = g[g >= 0]
        n = g.size
        if n == 0:
            continue
        obs = np.array([(g == 0).sum(), (g == 1).sum(), (g == 2).sum()],
                       dtype=float)
        p = (2 * obs[2] + obs[1]) / (2 * n)
        exp = np.array([(1 - p) ** 2, 2 * p * (1 - p), p ** 2]) * n
        keep = exp > 0
        if keep.sum() < 2:
            continue
        stat = float(((obs[keep] - exp[keep]) ** 2 / exp[keep]).sum())
        out[j] = chi2.sf(stat, 1)
    return out


def het_rate(G: np.ndarray) -> np.ndarray:
    """Per-sample heterozygosity over called genotypes."""
    called = (G >= 0).sum(axis=1)
    return (G == 1).sum(axis=1) / np.maximum(called, 1)


def _het_keep(G: np.ndarray, sd: float,
              ancestry: np.ndarray | None) -> np.ndarray:
    """Flag het-rate outliers as |z| > sd, computed within ancestry
    group when labels exist — a pooled z would mistake a real
    population's het level for contamination."""
    rate = het_rate(G)
    if ancestry is not None and len(ancestry) == len(rate) \
            and len(np.unique(ancestry[ancestry >= 0])) > 0:
        keep = np.ones(len(rate), bool)
        for a in np.unique(ancestry):
            if a < 0:
                continue
            m = ancestry == a
            mu, s = rate[m].mean(), rate[m].std()
            if s > 0:
                keep[m] = np.abs(rate[m] - mu) / s <= sd
        return keep
    mu, s = rate.mean(), rate.std()
    return np.abs(rate - mu) / max(s, 1e-12) <= sd


def apply_qc(G: np.ndarray, samples: np.ndarray, q: dict,
             ancestry: np.ndarray | None = None):
    """Apply filters in order; returns filtered G, samples, kept sample
    mask, kept variant indices (into the input), and the waterfall."""
    s_in, v_in = G.shape
    steps = []

    miss_s = (G < 0).mean(axis=1)
    keep_s = miss_s <= q["max_sample_missing"]
    steps.append({"name": "sample_missingness", "axis": "sample",
                  "removed": int(s_in - keep_s.sum()),
                  "kept": int(keep_s.sum())})
    G, samples = G[keep_s], samples[keep_s]
    if ancestry is not None:
        ancestry = ancestry[keep_s]

    if q.get("het_sd", 0) > 0:
        keep_h = _het_keep(G, q["het_sd"], ancestry)
        steps.append({"name": "heterozygosity", "axis": "sample",
                      "removed": int((~keep_h).sum()),
                      "kept": int(keep_h.sum())})
        G, samples = G[keep_h], samples[keep_h]
        if ancestry is not None:
            ancestry = ancestry[keep_h]
        # compose sample mask: keep_s gets het-rejects marked False
        idx = np.where(keep_s)[0]
        keep_s = keep_s.copy()
        keep_s[idx[~keep_h]] = False

    idx_v = np.arange(v_in)
    for name, mask in (
        ("variant_missingness",
         (G < 0).mean(axis=0) <= q["max_variant_missing"]),
        ("min_maf", None),
        ("hwe", None),
    ):
        if mask is None:
            if name == "min_maf":
                keep = maf(G) >= q["min_maf"]
            else:
                if q.get("hwe_p_min", 0) <= 0:
                    continue
                hp = hwe_p(G)
                keep = np.isnan(hp) | (hp >= q["hwe_p_min"])
        else:
            keep = mask
        steps.append({"name": name, "axis": "variant",
                      "removed": int((~keep).sum()),
                      "kept": int(keep.sum())})
        G, idx_v = G[:, keep], idx_v[keep]

    wf = {
        "samples_in": int(s_in), "samples_out": int(len(samples)),
        "variants_in": int(v_in), "variants_out": int(G.shape[1]),
        "steps": steps,
        "thresholds": {k: q[k] for k in
                       ("max_sample_missing", "max_variant_missing",
                        "min_maf", "hwe_p_min") if k in q},
    }
    if q.get("het_sd", 0) > 0:
        wf["thresholds"]["het_sd"] = q["het_sd"]
    assert wf["samples_in"] == sum(
        s["removed"] for s in wf["steps"] if s["axis"] == "sample"
    ) + wf["samples_out"]
    assert wf["variants_in"] == sum(
        s["removed"] for s in wf["steps"] if s["axis"] == "variant"
    ) + wf["variants_out"]
    return G, samples, keep_s, idx_v, wf


def main() -> None:
    cfg = load_config()
    data_dir = Path(cfg["paths"]["data_dir"])
    results_dir = Path(cfg["paths"]["results_dir"])
    z = np.load(data_dir / "raw" / "cohort.npz", allow_pickle=True)

    G, samples, keep_s, idx_v, wf = apply_qc(
        z["G"], z["samples"], cfg["qc"],
        ancestry=z["ancestry"] if "ancestry" in z else None)

    filt = dict(z)
    filt["G"] = G
    filt["samples"] = samples
    for key in ("vid", "chrom", "pos", "a1", "a2"):
        filt[key] = z[key][idx_v]
    for key in ("y", "cov", "ancestry"):
        if key in z:
            filt[key] = z[key][keep_s]

    out = data_dir / "processed"
    out.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out / "filtered.npz", **filt)
    results_dir.mkdir(parents=True, exist_ok=True)
    with open(results_dir / "qc_waterfall.json", "w") as f:
        json.dump(wf, f, indent=2)


if __name__ == "__main__":
    sys.exit(main())
