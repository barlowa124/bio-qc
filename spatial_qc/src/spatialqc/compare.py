"""Compare filtering strategies on one loaded sample.

For each strategy: spots retained, counts retained, genes still
detected, median depth of kept spots. When scanpy is importable, each
filtered matrix is Leiden-clustered and the report adds n_clusters and
modularity — a downstream-utility proxy, not a ground-truth score.
"""

from __future__ import annotations

import numpy as np

from . import filters, metrics


def compare(sample: dict, run_clusters: bool = True) -> dict:
    met = metrics.spot_metrics(sample)
    m = sample["matrix"]
    genes = np.asarray(sample["genes"])
    rows = {}
    for name, fn in filters.STRATEGIES.items():
        keep = fn(met)
        sub = m[keep]
        row = {
            "spots_kept": int(keep.sum()),
            "frac_spots": round(float(keep.mean()), 4),
            "counts_kept_frac": round(
                float(sub.sum() / m.sum()) if m.sum() else 0.0, 4),
            "genes_detected": int((np.asarray(sub.sum(axis=0))
                                   .ravel() > 0).sum()),
            "median_counts_kept":
                float(np.median(met["total_counts"][keep]))
                if keep.any() else 0.0,
            "median_mito_kept":
                float(np.median(met["mito_frac"][keep]))
                if keep.any() else 0.0,
        }
        if run_clusters and keep.sum() >= 50:
            row.update(_cluster(sub, genes))
        rows[name] = row
    return {"n_spots_total": int(m.shape[0]),
            "strategies": rows}


def _cluster(sub, genes) -> dict:
    try:
        import scanpy as sc
        import anndata
    except ImportError:
        return {"clusters": "skipped (scanpy unavailable)"}
    ad = anndata.AnnData(sub.astype(np.float32))
    ad.var_names = genes
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)
    sc.pp.scale(ad, max_value=10)
    sc.pp.neighbors(ad, n_neighbors=10)
    try:
        sc.tl.leiden(ad, resolution=0.5, key_added="cl",
                     flavor="igraph", n_iterations=2, directed=False)
    except TypeError:                    # scanpy<1.10
        sc.tl.leiden(ad, resolution=0.5, key_added="cl",
                     n_iterations=2, directed=False)
    mod = None
    try:
        import igraph
        import scipy.sparse
        coo = scipy.sparse.coo_matrix(ad.obsp["connectivities"])
        g = igraph.Graph()
        g.add_vertices(ad.n_obs)
        g.add_edges(list(zip(coo.row.tolist(), coo.col.tolist())))
        mod = g.modularity(ad.obs["cl"].cat.codes.tolist())
    except Exception:
        pass
    return {"n_clusters": int(ad.obs["cl"].nunique()),
            "leiden_modularity": (round(mod, 4)
                                  if isinstance(mod, float) else None)}
