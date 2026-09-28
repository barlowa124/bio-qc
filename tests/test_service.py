import io as stdio

import numpy as np
import pandas as pd
import pytest

from cytof_qc import service


def _two_pop_events(n=800, seed=0):
    rng = np.random.default_rng(seed)
    a = rng.normal(0, 1, (n, 3))
    b = rng.normal(6, 1, (n, 3))
    events = pd.DataFrame(
        np.vstack([a, b]), columns=["CD3", "CD19", "CD45"]
    )
    events["Time"] = np.arange(2 * n, dtype=float)  # QC channel, not clustered
    return events


def test_analyze_events_clusters_two_populations():
    report = service.analyze_events(_two_pop_events())
    assert report["n_events"] == 1600
    assert report["n_clusters_found"] >= 2
    assert report["phenotypic_channels"] == ["CD3", "CD19", "CD45"]
    assert report["qc_only_channels"] == ["Time"]
    assert len(report["embedding"]) == 1600
    assert len(report["embedding"][0]) == 2
    assert len(report["embedding_clusters"]) == 1600
    # Leiden over-partitions uniform blobs at default resolution, so the
    # check is separation: blob-A and blob-B events share almost no clusters
    labels = pd.Series(report["embedding_clusters"])
    a_clusters = set(labels.iloc[:800])
    b_clusters = set(labels.iloc[800:])
    overlap = a_clusters & b_clusters
    shared = sum(
        (labels == c).sum() for c in overlap
    )
    assert shared < 0.1 * len(labels)


def test_analyze_events_downsamples_embedding(monkeypatch):
    monkeypatch.setattr(service, "MAX_EMBED_POINTS", 100)
    report = service.analyze_events(_two_pop_events())
    assert len(report["embedding"]) == 100
    assert report["n_events"] == 1600


def test_analyze_events_rejects_oversize(monkeypatch):
    monkeypatch.setattr(service, "MAX_EVENTS", 10)
    with pytest.raises(ValueError, match="cap"):
        service.analyze_events(_two_pop_events())


def _client():
    from fastapi.testclient import TestClient

    return TestClient(service.create_app())


def test_health_endpoint():
    resp = _client().get("/health")
    assert resp.status_code == 200
    assert resp.json()["ok"] is True


def test_analyze_endpoint_accepts_csv():
    df = _two_pop_events()
    buf = stdio.BytesIO(df.to_csv(index=False).encode())
    resp = _client().post(
        "/api/analyze", files=[("files", ("events.csv", buf, "text/csv"))]
    )
    assert resp.status_code == 200, resp.text
    report = resp.json()
    assert report["n_events"] == 1600
    assert report["n_clusters_found"] >= 2
    assert report["acquisition_drift"]
    # uploads have no manual gates — no agreement claims
    assert "agreement" not in report


def test_analyze_endpoint_rejects_bad_file():
    resp = _client().post(
        "/api/analyze",
        files=[("files", ("x.bin", stdio.BytesIO(b"\x00\x01"), "app/octet"))],
    )
    assert resp.status_code == 400


def test_analyze_endpoint_rejects_empty_upload():
    resp = _client().post("/api/analyze", files=[])
    assert resp.status_code in (400, 422)
