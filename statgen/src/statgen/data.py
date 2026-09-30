"""Cohort stage: load or generate the raw genotype matrix.

Modes (config dataset.mode):
  demo - deterministic synthetic cohort with two ancestry groups and
         planted causal variants at recorded effect sizes. The planted
         truth is written next to the cohort so association recovery
         is auditable. Demo numbers are pipeline exercises, not
         findings.
  vcf  - biallelic SNPs from a local VCF (plain or gzipped, truncated
         streams tolerated). No public GWAS cohort ships a phenotype,
         so a trait is simulated on the real genotypes with a recorded
         seed; that is stated in the report and the truth file.

Matrix layout in cohort.npz:
  G        int8 S x V dosage of the counted allele; -1 = missing
  samples  sample ids
  vid/chrom/pos/a1/a2  per-variant fields
  y        phenotype (float)
  cov      S x K covariate matrix, cov_names column labels
  ancestry group labels when known (demo only)
"""

from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path

import numpy as np

from statgen.util import load_config


def _balding_nichols(rng, p: np.ndarray, fst: float) -> np.ndarray:
    """Draw a subpopulation allele frequency for each variant."""
    a = p * (1.0 - fst) / fst
    b = (1.0 - p) * (1.0 - fst) / fst
    return np.clip(rng.beta(a, b), 1e-4, 1.0 - 1e-4)


def generate_demo(cfg: dict) -> tuple[dict, dict]:
    d = cfg["dataset"]["demo"]
    rng = np.random.default_rng(d["seed"])
    n, v = d["n_samples"], d["n_variants"]
    share = d["ancestry_shares"]
    ancestry = rng.choice(len(share), size=n, p=share)

    p0 = rng.uniform(0.05, 0.45, size=v)
    p_anc = np.stack(
        [_balding_nichols(rng, p0, d["fst"]) for _ in share], axis=1)
    G = rng.binomial(2, p_anc[:, ancestry]).astype(np.int8).T  # S x V
    miss = rng.random(G.shape) < d["missing_rate"]
    G[miss] = -1

    cov = np.column_stack([
        rng.binomial(1, 0.5, size=n).astype(float),
        rng.normal(50.0, 10.0, size=n),
    ])
    cov_names = ["sex", "age"]

    k = d["n_causal"]
    causal_idx = np.sort(rng.choice(v, size=k, replace=False))
    betas = rng.uniform(d["beta_lo"], d["beta_hi"], size=k) * rng.choice(
        [-1.0, 1.0], size=k)
    Gf = np.where(G < 0, np.nan, G).astype(float)
    maf = np.nanmean(Gf, axis=0) / 2.0
    signal = np.nansum(Gf[:, causal_idx] * betas, axis=1)
    signal = np.nan_to_num(signal)
    anc_shift = d["ancestry_shift"] * ancestry.astype(float)
    var_g = np.var(signal) + np.var(anc_shift)
    noise_sd = np.sqrt(var_g * (1.0 - d["h2"]) / d["h2"])
    y = signal + anc_shift + rng.normal(0.0, noise_sd, size=n)

    vids = np.array([f"rs{100000 + i}" for i in range(v)])
    cohort = {
        "G": G, "samples": np.array([f"DG{i:04d}" for i in range(n)]),
        "vid": vids, "chrom": np.full(v, "synthetic"),
        "pos": np.arange(1, v + 1),
        "a1": np.full(v, "A"), "a2": np.full(v, "G"),
        "y": y.astype(float), "cov": cov,
        "cov_names": np.array(cov_names),
        "ancestry": ancestry.astype(np.int8),
    }
    truth = {
        "kind": "synthetic",
        "seed": d["seed"],
        "n_causal": k,
        "causal": [
            {"vid": str(vids[i]), "pos": int(i + 1), "beta": float(b),
             "maf": float(maf[i])}
            for i, b in zip(causal_idx, betas)
        ],
        "ancestry_shares": share, "fst": d["fst"],
        "ancestry_shift": d["ancestry_shift"], "h2_target": d["h2"],
        "noise_sd": float(noise_sd),
    }
    return cohort, truth


def _iter_vcf(path: str):
    """Yield (chrom, pos, vid, ref, alt, dosage array) for biallelic SNPs.

    Tolerates truncated gzip streams: a range-fetched prefix of a bgzip
    file ends mid-block and simply stops yielding.
    """
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", errors="replace") as f:
        samples: list[str] = []
        try:
            for line in f:
                if line.startswith("##"):
                    continue
                if line.startswith("#CHROM"):
                    samples = line.rstrip("\n").split("\t")[9:]
                    continue
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 10 or not samples:
                    continue
                ref, alt = parts[3], parts[4]
                if len(ref) != 1 or len(alt) != 1:
                    continue
                gts = parts[9:]
                dos = np.full(len(samples), -1, dtype=np.int8)
                for i, g in enumerate(gts):
                    gt = g.split(":", 1)[0].replace("|", "/")
                    if "." in gt:
                        continue
                    try:
                        dos[i] = sum(int(a) for a in gt.split("/"))
                    except ValueError:
                        continue
                yield parts[0], int(parts[1]), parts[2], ref, alt, dos
        except (EOFError, gzip.BadGzipFile, OSError):
            return


def vcf_samples(path: str) -> list[str]:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", errors="replace") as f:
        for line in f:
            if line.startswith("#CHROM"):
                return line.rstrip("\n").split("\t")[9:]
            if not line.startswith("#"):
                break
    return []


def load_vcf(cfg: dict) -> tuple[dict, dict]:
    d = cfg["dataset"]["vcf"]
    rng = np.random.default_rng(d["trait_seed"])
    samples = np.array(vcf_samples(d["path"]))
    rows = []
    for rec in _iter_vcf(d["path"]):
        chrom, pos, vid, ref, alt, dos = rec
        if vid in ("", "."):  # VCF "." = no rsID assigned
            vid = f"{chrom}:{pos}"
        rows.append((chrom, pos, vid, ref, alt, dos))
        if len(rows) >= d["max_variants"]:
            break
    if not rows:
        raise ValueError(f"no biallelic SNPs parsed from {d['path']}")

    n = len(rows[0][5])
    if len(samples) != n:
        samples = np.array([f"S{i}" for i in range(n)])
    G = np.stack([r[5] for r in rows], axis=1)
    meta = np.array(rows, dtype=object)

    # Optional sample panel (1000 Genomes integrated_call_samples):
    # superpopulation labels for stratification diagnostics and low-
    # differentiation causal-variant placement, plus sex as a real
    # covariate.
    ancestry = None
    cov = np.zeros((n, 0))
    cov_names = np.array([])
    panel = d.get("panel_path")
    if panel and Path(panel).exists():
        super_pop = {}
        sex = {}
        with open(panel) as f:
            for i, line in enumerate(f):
                if i == 0:
                    continue
                p = line.split("\t")
                if len(p) >= 4:
                    super_pop[p[0]] = p[2]
                    sex[p[0]] = p[3].strip()
        labels = sorted({super_pop.get(s) for s in samples} - {None})
        if labels:
            code = {lab: i for i, lab in enumerate(labels)}
            ancestry = np.array(
                [code.get(super_pop.get(s), -1) for s in samples],
                dtype=np.int8)
            cov = np.array(
                [1.0 if sex.get(s) == "male" else 0.0 for s in samples]
            ).reshape(-1, 1)
            cov_names = np.array(["sex_male"])

    # Simulated trait on real genotypes: a few planted variants + noise.
    # Causals come from common variants (pooled maf >= 0.05) so they can
    # survive QC, and — when a panel is present — from variants with
    # limited continental frequency spread so the trait stays
    # unconfounded (set dataset.vcf.ancestry_shift to deliberately
    # reintroduce stratification).
    Gm = np.where(G < 0, np.nan, G).astype(float)
    pooled_maf = np.nanmean(Gm, axis=0) / 2.0
    common = np.where(pooled_maf >= 0.05)[0]
    spread_ok = common
    if ancestry is not None and (ancestry >= 0).all():
        spread = np.full(G.shape[1], np.inf)
        for j in common:
            freqs = [np.nanmean(Gm[ancestry == a, j]) / 2.0
                     for a in np.unique(ancestry)]
            spread[j] = max(freqs) - min(freqs)
        spread_ok = common[spread[common] <= d.get("max_anc_spread", 0.15)]
    v = G.shape[1]
    k = min(d["n_causal"], len(spread_ok))
    causal_idx = np.sort(rng.choice(spread_ok, size=k, replace=False))
    betas = rng.uniform(d["beta_lo"], d["beta_hi"], size=k) * rng.choice(
        [-1.0, 1.0], size=k)
    Gf = np.where(G < 0, np.nan, G).astype(float)
    signal = np.nan_to_num(np.nansum(Gf[:, causal_idx] * betas, axis=1))
    shift = d.get("ancestry_shift", 0.0)
    if ancestry is not None and shift:
        signal = signal + shift * (
            ancestry.astype(float) - ancestry.astype(float).mean())
    noise_sd = np.sqrt(np.var(signal) * (1.0 - d["h2"]) / d["h2"])
    y = signal + rng.normal(0.0, noise_sd, size=n)

    cohort = {
        "G": G, "samples": samples,
        "vid": meta[:, 2].astype(str),
        "chrom": meta[:, 0].astype(str),
        "pos": meta[:, 1].astype(int),
        "a1": meta[:, 4].astype(str),
        "a2": meta[:, 3].astype(str),
        "y": y.astype(float),
        "cov": cov, "cov_names": cov_names,
    }
    if ancestry is not None:
        cohort["ancestry"] = ancestry
    truth = {
        "kind": "real_genotypes_simulated_trait",
        "source": d["path"], "trait_seed": d["trait_seed"],
        "ancestry_shift": shift,
        "n_causal": k,
        "causal": [
            {"vid": str(meta[i, 2]), "pos": int(meta[i, 1]),
             "beta": float(b)}
            for i, b in zip(causal_idx, betas)
        ],
        "h2_target": d["h2"], "noise_sd": float(noise_sd),
    }
    return cohort, truth


def main() -> None:
    cfg = load_config()
    out_dir = Path(cfg["paths"]["data_dir"]) / "raw"
    out_dir.mkdir(parents=True, exist_ok=True)

    mode = cfg["dataset"]["mode"]
    if mode == "demo":
        cohort, truth = generate_demo(cfg)
    elif mode == "vcf":
        cohort, truth = load_vcf(cfg)
    else:
        raise ValueError(f"unknown dataset.mode: {mode}")

    np.savez_compressed(out_dir / "cohort.npz", **cohort)
    Path(cfg["paths"]["results_dir"]).mkdir(parents=True, exist_ok=True)
    with open(cfg["paths"]["results_dir"] + "/planted_truth.json",
              "w") as f:
        json.dump(truth, f, indent=2)


if __name__ == "__main__":
    sys.exit(main())
