import numpy as np
import pandas as pd

from cytof_qc import mapping


def test_agreement_perfect_and_random():
    true = np.array(["a"] * 50 + ["b"] * 50)
    assert mapping.agreement(true, true)["adjusted_rand_index"] == 1.0
    rng = np.random.default_rng(0)
    rand = rng.permutation(true)
    assert mapping.agreement(true, rand)["adjusted_rand_index"] < 0.1


def test_match_clusters_swapped_labels():
    # clusters "0"/"1" recover populations "T"/"B" with names swapped
    true = pd.Series(["T"] * 40 + ["B"] * 40)
    pred = pd.Series(["0"] * 40 + ["1"] * 40)
    table = mapping.contingency(true, pred)
    m = mapping.match_clusters_to_populations(table)
    assert m == {"0": "T", "1": "B"}


def test_per_population_report_partial_split():
    # population "T" split across two clusters: recall < 1 for dominant
    true = pd.Series(["T"] * 10 + ["T"] * 5 + ["B"] * 20)
    pred = pd.Series(["0"] * 10 + ["1"] * 5 + ["1"] * 20)
    rep = mapping.per_population_report(true, pred)
    assert rep["T"]["recall"] == round(10 / 15, 4)
    assert rep["B"]["recall"] == 1.0
