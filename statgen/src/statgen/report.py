"""Report stage: summary.json, manhattan/qq figures, report.md, and the
claim check binding every number in the markdown to a recorded result
JSON (vendored claims.py — shared with scrna_qc, parity-pinned).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy

from statgen.claims import flatten_results, verify_markdown
from statgen.util import git_sha, load_config


def manhattan(df: pd.DataFrame, bonf: float, causal_ids: set | None,
              path: Path) -> None:
    fig, ax = plt.subplots(figsize=(9, 3.6))
    lp = -np.log10(df["p"].clip(lower=1e-300))
    ax.scatter(range(len(df)), lp, s=6, c="#4c72b0", linewidths=0)
    if causal_ids:
        m = df["vid"].isin(causal_ids)
        ax.scatter(np.where(m)[0], lp[m], s=14, c="#c44e52",
                   linewidths=0, label="planted causal", zorder=3)
        ax.legend(frameon=False, markerscale=1.6, fontsize=8)
    ax.axhline(-np.log10(bonf), color="#dd8452", ls="--", lw=1,
               label=None)
    ax.text(len(df) * 0.99, -np.log10(bonf), "Bonferroni",
            ha="right", va="bottom", fontsize=7, color="#dd8452")
    ax.set_ylabel(r"$-\log_{10}(p)$")
    ax.set_xlabel("variant")
    ax.set_xticks([])
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def qq(df: pd.DataFrame, lam: float, path: Path) -> None:
    obs = -np.log10(np.sort(df["p"].to_numpy()))
    n = len(obs)
    exp = -np.log10((np.arange(1, n + 1) - 0.5) / n)
    fig, ax = plt.subplots(figsize=(3.6, 3.6))
    ax.scatter(exp, obs, s=6, c="#4c72b0", linewidths=0)
    lim = max(exp[-1], obs[-1]) * 1.05
    ax.plot([0, lim], [0, lim], c="#c44e52", lw=1)
    ax.set_xlabel(r"expected $-\log_{10}(p)$")
    ax.set_ylabel(r"observed $-\log_{10}(p)$")
    ax.set_title(f"$\\lambda_{{GC}}$ = {lam:.3f}", fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def main() -> None:
    cfg = load_config()
    results_dir = Path(cfg["paths"]["results_dir"])

    wf = json.loads((results_dir / "qc_waterfall.json").read_text())
    strat = json.loads((results_dir / "strat_summary.json").read_text())
    assoc = json.loads((results_dir / "assoc_summary.json").read_text())
    truth = json.loads((results_dir / "planted_truth.json").read_text())
    df = pd.read_csv(results_dir / "assoc.tsv", sep="\t")

    causal_ids = {c["vid"] for c in truth.get("causal", [])}
    causal_pos = [c["pos"] for c in truth.get("causal", [])]
    causal_df = df[df["vid"].isin(causal_ids)]
    truth_beta = {c["vid"]: c["beta"] for c in truth.get("causal", [])}
    recovered = sum(
        1 for _, r in causal_df.iterrows()
        if r["p"] < assoc["bonferroni_alpha"])

    # LD accounting: non-causal hits within LD_WINDOW of a planted
    # causal are proxies of its signal, not independent discoveries.
    # Skipped on the demo cohort where positions are synthetic indices.
    hits = df[df["p"] < assoc["bonferroni_alpha"]]
    real_coords = cfg["dataset"]["mode"] != "demo"
    if real_coords:
        LD_WINDOW = 250_000
        n_hits_other = int(sum(
            1 for p in hits["pos"]
            if all(abs(p - cp) > LD_WINDOW for cp in causal_pos)))
    else:
        n_hits_other = int(len(hits) - recovered)
    beta_corr = None
    if len(causal_df) >= 2:
        beta_corr = float(np.corrcoef(
            causal_df["beta"],
            [truth_beta[v] for v in causal_df["vid"]])[0, 1])

    summary = {
        "pipeline": "statgen",
        "mode": cfg["dataset"]["mode"],
        "trait_kind": truth["kind"],
        "samples_in": wf["samples_in"], "samples_out": wf["samples_out"],
        "variants_in": wf["variants_in"],
        "variants_out": wf["variants_out"],
        "assoc_kind": assoc["kind"],
        "lambda_gc": assoc["lambda_gc"],
        "lambda_gc_no_cov": assoc["lambda_gc_no_cov"],
        "n_bonferroni_hits": assoc["n_bonferroni_hits"],
        "bonferroni_alpha": assoc["bonferroni_alpha"],
        "lead_hit": assoc["lead_hit"],
        "n_causal_planted": truth["n_causal"],
        "n_causal_post_qc": int(len(causal_df)),
        "n_causal_recovered": recovered,
        "causal_beta_corr": beta_corr,
        "n_hits_causal": int(recovered),
        "n_hits_ld_proxy": int(
            assoc["n_bonferroni_hits"] - recovered - n_hits_other),
        "n_hits_other": n_hits_other,
        "ld_window_kb": 250 if real_coords else 0,
        "pc1_ancestry_r": strat.get("pc1_ancestry_r"),
        "kinship_offdiag_max": strat.get("kinship_offdiag_max"),
        "n_related_pairs": strat.get("n_related_pairs"),
        "provenance": {
            "git": git_sha(),
            "python": sys.version.split()[0],
            "numpy": np.__version__, "scipy": scipy.__version__,
            "pandas": pd.__version__,
        },
    }
    with open(results_dir / "summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    manhattan(df, assoc["bonferroni_alpha"], causal_ids or None,
              results_dir / "manhattan.png")
    qq(df, assoc["lambda_gc"], results_dir / "qq.png")

    lead = assoc["lead_hit"]
    lines = [
        "# statgen report",
        "",
        f"Cohort mode `{summary['mode']}`; trait `{truth['kind']}`. "
        "Synthetic and simulated-trait numbers exercise the pipeline "
        "and are not genetic findings.",
        "",
        "## QC waterfall",
        "",
        f"{wf['samples_in']} -> {wf['samples_out']} samples; "
        f"{wf['variants_in']} -> {wf['variants_out']} variants.",
        "",
        "| filter | axis | removed | kept |",
        "|---|---|---:|---:|",
    ]
    for s in wf["steps"]:
        lines.append(
            f"| {s['name']} | {s['axis']} | {s['removed']} "
            f"| {s['kept']} |")
    ve = strat["var_explained"][0] * 100
    lines += [
        "",
        "## Stratification",
        "",
        f"{strat['n_pcs']} PCs computed on "
        f"{strat['n_variants_used']} variants; PC1 explains "
        f"{ve:.1f}% of genotype variance",
    ]
    if strat.get("n_related_pairs") is not None:
        lines += [
            "",
            "Relatedness scan (GRM off-diagonal after projecting out "
            f"the {strat['n_pcs']} PCs, so ancestry sharing does not "
            f"count as relatedness): max {strat['kinship_offdiag_max']:.3f}, "
            f"sd {strat['kinship_offdiag_sd']:.3f}, "
            f"{strat['n_related_pairs']} pairs above the "
            f"{strat['kinship_flag']} flag"
            + (" — with few regional variants the estimate is "
               "overdispersed; a high count reflects LD, not cryptic "
               "relatedness." if strat["n_variants_used"] < 50000
               and strat["n_related_pairs"] > 0 else "."),
        ]
    if strat.get("pc1_ancestry_r") is not None:
        lines[-1] += (
            f" and correlates with the ancestry label at "
            f"|r| = {abs(strat['pc1_ancestry_r']):.2f}")
    elif strat.get("pc1_ancestry_eta2") is not None:
        lines[-1] += (
            f"; {strat['n_ancestry_groups']} ancestry groups capture "
            f"{strat['pc1_ancestry_eta2'] * 100:.1f}% of its variance "
            f"(eta^2)")
    lines[-1] += "."
    lines += [
        "",
        "## Association",
        "",
        f"{assoc['kind']} test, {assoc['n_variants']} variants on "
        f"{assoc['n_samples']} samples with "
        f"{assoc['n_covariates']} covariates "
        f"({assoc['n_pcs_as_cov']} PCs).",
        "",
        f"Genomic control $\\lambda$ = {assoc['lambda_gc']:.3f} "
        f"with covariates vs {assoc['lambda_gc_no_cov']:.3f} without; "
        f"{assoc['n_bonferroni_hits']} variants pass the Bonferroni "
        f"threshold {assoc['bonferroni_alpha']:.2e}.",
        "",
        "$\\lambda$ by covariate PCs: "
        + ", ".join(f"{k.split('_')[1]} PCs -> {v:.3f}"
                    for k, v in sorted(
                        assoc["lambda_by_pcs"].items(),
                        key=lambda kv: int(kv[0].split("_")[1])))
        + ". Effect sizes are per copy of the counted allele (a1).",
        "",
        f"Lead hit {lead['vid']} (pos {lead['pos']}): "
        f"beta {lead['beta']:.3f} (se {lead['se']:.3f}), "
        f"p = {lead['p']:.2e}.",
        "",
        f"Planted causal variants: {summary['n_causal_recovered']} "
        f"of {summary['n_causal_post_qc']} surviving QC recovered at "
        f"Bonferroni ({summary['n_causal_planted']} planted"
        + (f"; {summary['n_causal_planted'] - summary['n_causal_post_qc']}"
           " failed pooled HWE — the Wahlund effect of mixing groups"
           if summary["n_causal_planted"] != summary["n_causal_post_qc"]
           else "")
        + (f"; beta correlation {beta_corr:.2f}" if beta_corr
           is not None else "")
        + ".",
        "",
        f"Of the {assoc['n_bonferroni_hits']} Bonferroni hits, "
        + (f"{summary['n_hits_ld_proxy']} sit within 250 kb of a "
           f"planted causal (LD proxies) and " if real_coords else "")
        + f"{summary['n_hits_other']} are "
        + ("unexplained — marginal tests cannot separate correlated "
           "variants, so per-variant attribution needs fine-mapping."
           if real_coords else "non-causal."),
        "",
        "![manhattan](manhattan.png)",
        "",
        "![qq](qq.png)",
        "",
        "## Limits",
        "",
        "Single-marker tests on a small cohort; no LD-aware "
        "fine-mapping or imputation. Simulated traits use recorded "
        "seeds and planted effects only to exercise recovery.",
    ]
    md = "\n".join(lines) + "\n"
    (results_dir / "report.md").write_text(md)

    pool: dict = {}
    for name in ("qc_waterfall", "strat_summary", "assoc_summary",
                 "summary", "planted_truth"):
        src = json.loads((results_dir / f"{name}.json").read_text())
        for k, v in flatten_results(src).items():
            pool[f"{name}.{k}"] = v
    verdict = verify_markdown(md, pool)
    verdict["verified_by"] = "statgen.claims (vendored, parity-pinned)"
    with open(results_dir / "claims_check.json", "w") as f:
        json.dump(verdict, f, indent=2)
    if not verdict["passed"]:
        print("UNBOUND CLAIMS:",
              [u["token"] for u in verdict["unbound_claims"]])


if __name__ == "__main__":
    sys.exit(main())
