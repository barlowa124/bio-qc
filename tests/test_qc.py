import numpy as np
import pandas as pd

from cytof_qc import qc


def _events(n=100):
    rng = np.random.default_rng(0)
    return pd.DataFrame(
        {
            "CD45": rng.uniform(1, 100, n),
            "CD3": rng.uniform(0, 50, n),
            "population": np.repeat(["T", "B"], n // 2),
        }
    )


def test_channel_qc_flags_injected_negatives():
    ev = _events()
    ev.loc[ev.index[:50], "CD3"] = -1.0  # 50% negative in CD3
    report = qc.channel_qc(ev, ["CD45", "CD3"])
    assert report["CD3"]["flag_negative"] is True
    assert report["CD45"]["flag_negative"] is False
    assert report["CD3"]["negative_fraction"] == 0.5


def test_channel_qc_flags_dead_channel():
    ev = _events()
    ev["CD3"] = 0.0
    report = qc.channel_qc(ev, ["CD3"])
    assert report["CD3"]["flag_zero"] is True


def test_population_qc_flags_sparse():
    ev = _events()
    ev.loc[ev.index[:3], "population"] = "rare"
    report = qc.population_qc(ev, min_events=20)
    assert report["rare"]["flag_sparse"] is True
    assert report["T"]["flag_sparse"] is False
