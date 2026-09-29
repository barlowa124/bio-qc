import io as stdio
import time

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
    client = _client()
    resp = client.post(
        "/api/analyze", files=[("files", ("events.csv", buf, "text/csv"))]
    )
    assert resp.status_code == 200, resp.text
    job = resp.json()
    assert job["status"] == "running"
    assert job["n_events"] == 1600

    report = None
    for _ in range(120):
        jr = client.get(f"/api/jobs/{job['job_id']}")
        body = jr.json()
        if body.get("status") == "running":
            import time

            time.sleep(0.5)
            continue
        report = body
        break
    assert report is not None, "job never finished"
    assert report["n_events"] == 1600
    assert report["n_clusters_found"] >= 2
    assert report["acquisition_drift"]
    assert report["timing"]["analyze_ms"] > 0
    # uploads have no manual gates — no agreement claims
    assert "agreement" not in report


def test_unknown_job_404():
    resp = _client().get("/api/jobs/deadbeefdead")
    assert resp.status_code == 404


def test_hung_job_expires_to_error():
    service._JOBS["hungjob12345"] = {
        "status": "running",
        "started": time.time() - service.MAX_JOB_RUNTIME_S - 1,
        "n_events": 100,
    }
    try:
        resp = _client().get("/api/jobs/hungjob12345")
        assert resp.status_code == 200
        body = resp.json()
        assert body["status"] == "error"
        assert "runtime" in body["error"]
    finally:
        service._JOBS.pop("hungjob12345", None)


def test_analyze_endpoint_rejects_tiny_upload():
    df = _two_pop_events(n=5)  # 10 events < MIN_EVENTS
    buf = stdio.BytesIO(df.to_csv(index=False).encode())
    resp = _client().post(
        "/api/analyze", files=[("files", ("tiny.csv", buf, "text/csv"))]
    )
    assert resp.status_code == 400
    assert "too few" in resp.json()["error"]


def test_analyze_events_rejects_tiny_input():
    with pytest.raises(ValueError, match="too few"):
        service.analyze_events(_two_pop_events(n=5))


def test_analyze_endpoint_rejects_bad_file():
    resp = _client().post(
        "/api/analyze",
        files=[("files", ("x.bin", stdio.BytesIO(b"\x00\x01"), "app/octet"))],
    )
    assert resp.status_code == 400


def test_analyze_endpoint_rejects_empty_upload():
    resp = _client().post("/api/analyze", files=[])
    assert resp.status_code in (400, 422)


def test_upload_over_file_cap_rejected(monkeypatch):
    # even a small file is over the cap once MAX_FILE_BYTES is lowered —
    # exercises the bounded-read path without needing a 200MB fixture
    monkeypatch.setattr(service, "MAX_FILE_BYTES", 64)
    df = _two_pop_events(n=40)
    buf = stdio.BytesIO(df.to_csv(index=False).encode())
    resp = _client().post(
        "/api/analyze", files=[("files", ("big.csv", buf, "text/csv"))]
    )
    assert resp.status_code == 400
    assert "file cap" in resp.json()["error"]


def test_busy_returns_429(monkeypatch):
    monkeypatch.setattr(service, "MAX_RUNNING_JOBS", 0)
    df = _two_pop_events(n=40)
    buf = stdio.BytesIO(df.to_csv(index=False).encode())
    resp = _client().post(
        "/api/analyze", files=[("files", ("ok.csv", buf, "text/csv"))]
    )
    assert resp.status_code == 429


def test_finished_jobs_expire_past_ttl():
    service._JOBS["oldjob000000"] = {
        "status": "done",
        "started": time.time() - service.JOB_TTL_S - 10,
        "n_events": 50,
        "report": {},
    }
    df = _two_pop_events(n=40)
    buf = stdio.BytesIO(df.to_csv(index=False).encode())
    resp = _client().post(
        "/api/analyze", files=[("files", ("ok.csv", buf, "text/csv"))]
    )
    assert resp.status_code == 200
    assert "oldjob000000" not in service._JOBS


def test_stats_endpoint_reports_log_summary(tmp_path, monkeypatch):
    log = tmp_path / "requests.jsonl"
    log.write_text(
        '\n'.join([
            '{"ts":1,"path":"/api/analyze","status":200,"ms":100,"n_events":500}',
            '{"ts":2,"path":"/api/analyze","status":400,"ms":20}',
            '{"ts":3,"path":"/health","status":200,"ms":2}',
        ]) + "\n"
    )
    monkeypatch.setattr(service, "_REQUEST_LOG", log)
    resp = _client().get("/api/stats")
    assert resp.status_code == 200
    body = resp.json()
    assert body["requests"] == 3
    assert body["by_status"]["200"] == 2
    assert body["by_status"]["400"] == 1
    assert body["latency_ms"]["p50"] == 20
    assert body["n_events"]["max"] == 500
