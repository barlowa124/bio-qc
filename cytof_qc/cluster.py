"""Unsupervised clustering on transformed events via anndata/scanpy."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc


def cluster_events(
    events_t: pd.DataFrame,
    channels: list[str],
    *,
    n_neighbors: int = 15,
    resolution: float = 1.0,
    random_state: int = 0,
) -> ad.AnnData:
    """Leiden-cluster events on the selected marker channels.

    Input is already arcsinh-transformed. Returns the AnnData with
    ``.obs["cluster"]``, UMAP coordinates in ``.obsm["X_umap"]``, and the
    original transformed values in ``.X``.
    """
    adata = ad.AnnData(events_t[channels].to_numpy(dtype=np.float32))
    adata.var_names = channels
    adata.obs_names = [str(i) for i in events_t.index]
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, random_state=random_state)
    sc.tl.leiden(
        adata,
        resolution=resolution,
        random_state=random_state,
        flavor="igraph",
        n_iterations=2,
        directed=False,
        key_added="cluster",
    )
    sc.tl.umap(adata, random_state=random_state)
    return adata
