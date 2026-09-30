"""Per-cell QC metrics: library size, detected genes, mitochondrial %.

Computed directly from X rather than delegated to a library helper so the
definitions are inspectable and testable against hand-worked fixtures.
"""

from __future__ import annotations

import numpy as np
import anndata as ad


def _row_sums(X) -> np.ndarray:
    return np.asarray(X.sum(axis=1)).ravel()


def _row_nnz(X) -> np.ndarray:
    return np.asarray((X > 0).sum(axis=1)).ravel()


def add_qc_metrics(adata: ad.AnnData, mt_prefix: str) -> ad.AnnData:
    X = adata.X
    total = _row_sums(X)
    adata.obs["total_counts"] = total
    adata.obs["n_genes_by_counts"] = _row_nnz(X)

    mt_mask = adata.var_names.str.startswith(mt_prefix)
    if mt_mask.any():
        mt_counts = _row_sums(X[:, mt_mask])
        adata.obs["pct_mt"] = np.where(total > 0, 100.0 * mt_counts / total, 0.0)
    else:
        adata.obs["pct_mt"] = np.zeros(adata.n_obs)
    return adata
