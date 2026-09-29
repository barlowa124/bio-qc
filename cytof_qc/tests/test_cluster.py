import numpy as np
import pandas as pd
import pytest

from cytof_qc import cluster


def test_cluster_events_recovers_two_blobs():
    rng = np.random.default_rng(0)
    a = rng.normal(2.0, 0.3, (300, 4))
    b = rng.normal(5.0, 0.3, (300, 4))
    events = pd.DataFrame(np.vstack([a, b]), columns=list("wxyz"))
    try:
        adata = cluster.cluster_events(
            events, list("wxyz"), n_neighbors=15, resolution=0.2
        )
    except ImportError:
        pytest.skip("igraph/leidenalg not installed")
    # clusters may overfragment at higher resolution, but every cluster
    # must be dominated by one blob
    obs = adata.obs["cluster"].to_numpy()
    assert len(np.unique(obs)) >= 2
    for cl in np.unique(obs):
        idx = np.where(obs == cl)[0]
        frac_first_blob = np.mean(idx < 300)
        assert frac_first_blob < 0.1 or frac_first_blob > 0.9
