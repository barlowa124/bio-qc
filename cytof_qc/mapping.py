"""Score unsupervised clusters against gated reference populations.

Reports agreement metrics (ARI, NMI) and a per-population view of how the
clusters recover each manual gate, without claiming either is ground
truth in an absolute sense.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score


def agreement(labels_true: np.ndarray, labels_pred: np.ndarray) -> dict:
    """ARI and NMI between two labelings of the same events."""
    return {
        "adjusted_rand_index": round(float(adjusted_rand_score(labels_true, labels_pred)), 4),
        "normalized_mutual_info": round(
            float(normalized_mutual_info_score(labels_true, labels_pred)), 4
        ),
    }


def contingency(labels_true: pd.Series, labels_pred: pd.Series) -> pd.DataFrame:
    """population x cluster count table."""
    return pd.crosstab(labels_true, labels_pred)


def match_clusters_to_populations(table: pd.DataFrame) -> dict[str, str]:
    """Map each cluster to its best population via the Hungarian algorithm.

    `table` is populations (rows) x clusters (cols); matching maximizes
    total shared events. Returns cluster -> population.
    """
    cost = -table.to_numpy(dtype=float)  # maximize overlap
    row_idx, col_idx = linear_sum_assignment(cost)
    clusters = table.columns.to_numpy()
    populations = table.index.to_numpy()
    return {str(clusters[c]): str(populations[r]) for r, c in zip(row_idx, col_idx)}


def per_population_report(
    labels_true: pd.Series, labels_pred: pd.Series
) -> dict[str, dict]:
    """For each gated population: dominant cluster, recall, and purity."""
    table = contingency(labels_true, labels_pred)
    report = {}
    for pop, row in table.iterrows():
        total = int(row.sum())
        if total == 0:
            continue
        dom_cluster = str(row.idxmax())
        recovered = int(row.max())
        col = table[dom_cluster]
        report[str(pop)] = {
            "n_events": total,
            "dominant_cluster": dom_cluster,
            "recall": round(recovered / total, 4),
            "purity": round(recovered / int(col.sum()), 4),
        }
    return report
