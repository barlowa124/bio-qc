"""Fidelity scoring: organoid clusters vs reference cell-type centroids.

For each organoid Leiden cluster:
  - mean log-normalized HVG profile -> Pearson correlation against every
    reference cell-type centroid -> best match (or "unmapped" below
    config.score.mapping_correlation_min)
  - detection of the matched type's marker panel

Aggregate fidelity signals (all reported per-population AND pooled, because
pooled scores hide missing/aberrant populations):
  - celltype_coverage: fraction of reference types hit by >=1 mapped cluster
  - unmapped_fraction: fraction of organoid cells in unmapped clusters
  - marker_detection: per type, fraction of expected markers detected in
    the best-matching cluster
  - qc_flags: pass/warn decisions vs config.score thresholds
"""

from __future__ import annotations

import json
import sys

import anndata as ad
import numpy as np
import pandas as pd
from organoid_qc.util import load_config
from scipy import sparse


def _mean_profile(adata: ad.AnnData, mask, genes) -> np.ndarray:
    X = adata.layers["lognorm"][np.asarray(mask)]
    if sparse.issparse(X):
        X = X.toarray()
    idx = adata.var_names.get_indexer(genes)
    missing = np.asarray(genes)[idx < 0]
    if len(missing):
        raise ValueError(
            f"{len(missing)} genes absent from var_names "
            f"(e.g. {list(missing[:3])}), run preprocess to align gene spaces"
        )
    return np.asarray(X[:, idx]).mean(axis=0)


def _corr(a: np.ndarray, b: np.ndarray) -> float:
    if a.std() == 0 or b.std() == 0:
        return 0.0
    return float(np.corrcoef(a, b)[0, 1])


def score(organoid: ad.AnnData, reference: ad.AnnData, cfg: dict) -> dict:
    ct_col = cfg["dataset"]["obs_columns"]["cell_type"]
    cl_col = cfg["dataset"]["obs_columns"]["cluster"]
    thresh = cfg["score"]

    hvg = organoid.var_names[organoid.var["highly_variable"]]
    ref_types = sorted(reference.obs[ct_col].unique())
    centroids = {
        ct: _mean_profile(reference, reference.obs[ct_col] == ct, hvg)
        for ct in ref_types
    }

    clusters = sorted(organoid.obs[cl_col].unique())
    cluster_rows = []
    for cl in clusters:
        mask = organoid.obs[cl_col] == cl
        prof = _mean_profile(organoid, mask, hvg)
        corrs = {ct: _corr(prof, c) for ct, c in centroids.items()}
        best_ct = max(corrs, key=corrs.get)
        best_corr = corrs[best_ct]
        mapped = best_corr >= thresh["mapping_correlation_min"]
        cluster_rows.append(
            {
                "cluster": cl,
                "n_cells": int(mask.sum()),
                "best_match": best_ct if mapped else None,
                "best_correlation": round(best_corr, 4),
                "mapped": bool(mapped),
                "correlations": {k: round(v, 4) for k, v in corrs.items()},
            }
        )

    n_cells = organoid.n_obs
    unmapped_cells = sum(
        r["n_cells"] for r in cluster_rows if not r["mapped"]
    )
    covered = {r["best_match"] for r in cluster_rows if r["mapped"]}
    coverage = len(covered) / len(ref_types)

    # Marker detection: expected markers reaching >= marker_expression_ratio
    # of their reference-type mean in the cluster(s) mapped to that type.
    # Uncovered types score 0.0 (no mapped cluster to detect them in).
    ratio = thresh["marker_expression_ratio"]
    marker_sets = organoid.uns.get("markers") or cfg["markers"]
    marker_detection = {}
    for ct, genes in marker_sets.items():
        genes = [
            g
            for g in genes
            if g in organoid.var_names and g in reference.var_names
        ]
        if not genes or ct not in centroids:
            marker_detection[ct] = None
            continue
        if ct not in covered:
            marker_detection[ct] = 0.0
            continue
        best = max(
            (r for r in cluster_rows if r["best_match"] == ct),
            key=lambda r: r["best_correlation"],
        )
        mask = organoid.obs[cl_col] == best["cluster"]
        org_mean = _mean_profile(organoid, mask, genes)
        ref_mean = _mean_profile(
            reference, reference.obs[ct_col] == ct, genes
        )
        detected = org_mean >= ratio * ref_mean
        marker_detection[ct] = round(float(detected.mean()), 4)

    detected_vals = [v for v in marker_detection.values() if v is not None]
    flags = {
        "celltype_coverage": bool(
            coverage >= thresh["min_celltype_coverage"]
        ),
        "unmapped_fraction": bool(
            unmapped_cells / n_cells <= thresh["max_unmapped_fraction"]
        ),
        # No assessable markers => fail, not vacuous pass: an empty panel
        # means the QC question could not be asked at all.
        "marker_detection": bool(
            detected_vals
            and np.mean(detected_vals) >= thresh["min_marker_detection"]
        ),
    }
    return {
        "n_organoid_cells": int(n_cells),
        "n_reference_cells": int(reference.n_obs),
        "reference_types": ref_types,
        "clusters": cluster_rows,
        "aggregate": {
            "celltype_coverage": round(coverage, 4),
            "unmapped_fraction": round(unmapped_cells / n_cells, 4),
            "mean_marker_detection": (
                round(float(np.mean(detected_vals)), 4)
                if detected_vals
                else None
            ),
            "per_type_marker_detection": marker_detection,
        },
        "qc_flags": flags,
        "qc_pass": all(flags.values()),
        "thresholds": thresh,
    }


def main() -> None:
    org_path, ref_path, out_path = sys.argv[1:4]
    cfg = load_config()
    result = score(ad.read_h5ad(org_path), ad.read_h5ad(ref_path), cfg)
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2)
    print(
        f"score: coverage {result['aggregate']['celltype_coverage']}, "
        f"unmapped {result['aggregate']['unmapped_fraction']}, "
        f"qc_pass={result['qc_pass']} -> {out_path}"
    )


if __name__ == "__main__":
    main()
