"""HTTP surface for the QC pipeline — upload events, get a report.

``analyze_events`` is the product core: any caller that can produce an
event table (uploaded FCS files, CSV exports, FlowJo workspaces) gets
channel QC, a drift gate, Leiden clusters, and a UMAP embedding back.
Uploaded data carries no manual gates, so the report is descriptive —
cluster structure plus QC flags — not agreement scoring.
"""

import json
import logging
import os
import tempfile
import threading
import time
import uuid
from collections import deque
from pathlib import Path

import numpy as np
import pandas as pd

from . import io, qc, transform

MAX_UPLOAD_FILES = 200
MAX_EVENTS = 500_000
MAX_EMBED_POINTS = 20_000
MAX_FILE_BYTES = 200_000_000  # per uploaded file, before parse
MIN_EVENTS = 16               # clustering needs >= n_neighbors events

# State lives under CYTOF_STATE_DIR, a mounted Railway volume path, or
# the repo dir — in that order. On a volume it survives redeploys;
# without one it degrades to per-container state as before.
_STATE_DIR = Path(
    os.environ.get("CYTOF_STATE_DIR")
    or os.environ.get("RAILWAY_VOLUME_MOUNT_PATH")
    or Path(__file__).resolve().parent.parent
)
_REQUEST_LOG = _STATE_DIR / "service_requests.jsonl"
_JOBS_FILE = _STATE_DIR / "jobs.jsonl"

# analysis jobs: clustering a 50k-event upload takes ~90s on a shared
# CPU — too long to hold an HTTP request open. POST returns a job id,
# GET /api/jobs/{id} polls it. Job transitions append to _JOBS_FILE so
# a restart can replay them; a job 'running' at load died with the
# process and is marked error.
_JOBS = {}
_JOBS_LOCK = threading.Lock()
MAX_RUNNING_JOBS = 2
JOB_TTL_S = 3600
MAX_JOB_RUNTIME_S = 900

# per-IP submit throttle — the endpoint is open, so without a cap a
# single client can hold the CPU busy indefinitely
_RL_LOCK = threading.Lock()
_SUBMIT_TIMES = {}  # client ip -> deque of recent POST timestamps
RATE_LIMIT_WINDOW_S = 3600
RATE_LIMIT_MAX = 12


def _rate_limit_ok(client_ip):
    """True while client_ip is under RATE_LIMIT_MAX submits per window."""
    now = time.time()
    with _RL_LOCK:
        times = _SUBMIT_TIMES.setdefault(client_ip, deque())
        while times and now - times[0] > RATE_LIMIT_WINDOW_S:
            times.popleft()
        if len(times) >= RATE_LIMIT_MAX:
            return False
        times.append(now)
        return True


def _persist_job(job_id, job):
    """Append the job's current state to the durable log."""
    try:
        with _JOBS_FILE.open("a") as fh:
            fh.write(json.dumps({"job_id": job_id, **job}, default=str) + "\n")
    except OSError:
        logging.getLogger("cytof_qc.service").warning("job log write failed")


def _load_jobs():
    """Replay the durable job log into _JOBS, then compact and rewrite.
    A job still 'running' at load died with the process — mark it error
    rather than resuming silently."""
    if not _JOBS_FILE.exists():
        return
    try:
        for line in _JOBS_FILE.read_text().splitlines():
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                continue
            jid = rec.pop("job_id", None)
            if jid:
                _JOBS[jid] = rec
    except OSError:
        return
    for j in _JOBS.values():
        if j["status"] == "running":
            j.update(status="error", error="interrupted by restart")
    _compact_jobs()


def _compact_jobs(now=None):
    """Drop finished jobs past TTL, keep last 50 finished + all running;
    rewrite the durable log so pruned entries don't resurrect."""
    now = time.time() if now is None else now
    running = {k: v for k, v in _JOBS.items() if v["status"] == "running"}
    finished = sorted(
        (kv for kv in _JOBS.items() if kv[1]["status"] != "running"),
        key=lambda kv: kv[1]["started"],
    )
    finished = [
        kv for kv in finished if now - kv[1]["started"] <= JOB_TTL_S
    ][-50:]
    _JOBS.clear()
    _JOBS.update(finished)
    _JOBS.update(running)
    try:
        with _JOBS_FILE.open("w") as fh:
            for jid, j in _JOBS.items():
                fh.write(json.dumps({"job_id": jid, **j}, default=str) + "\n")
    except OSError:
        pass


def _expire_hung_locked(now=None):
    """Mark jobs still 'running' past the runtime cap as failed. Call
    with _JOBS_LOCK held. Guards against a hung worker permanently
    holding a MAX_RUNNING_JOBS slot — the daemon thread may still be
    stuck, but slot bookkeeping no longer waits on it."""
    now = time.time() if now is None else now
    for jid, j in _JOBS.items():
        if j["status"] == "running" and now - j["started"] > MAX_JOB_RUNTIME_S:
            j.update(status="error", error="job exceeded runtime cap")
            _persist_job(jid, j)


def analyze_events(events: pd.DataFrame) -> dict:
    """Run QC + clustering on an ungated event table.

    Every numeric column is a candidate channel; known acquisition/QC
    channels are split out and never drive clustering, the same rule the
    benchmark pipeline applies.
    """
    if len(events) > MAX_EVENTS:
        msg = f"{len(events)} events exceeds the {MAX_EVENTS} cap"
        raise ValueError(msg)
    numeric = events.select_dtypes(include="number")
    all_channels = list(numeric.columns)
    channels, qc_only = qc.split_channels(all_channels)
    if not channels:
        msg = "no phenotypic channels found in upload"
        raise ValueError(msg)
    if len(events) < MIN_EVENTS:
        msg = f"{len(events)} events is too few to cluster (need >= {MIN_EVENTS})"
        raise ValueError(msg)

    transformed = transform.arcsinh(numeric)
    channel_report = qc.channel_qc(
        events.assign(**{c: numeric[c] for c in numeric.columns}),
        all_channels,
    )
    drift_report = qc.acquisition_qc(transformed, all_channels)

    from . import cluster

    t_cluster = time.time()
    # cluster on all events; the embedding is display-only, so it is
    # computed on the subsample the browser actually receives
    adata = cluster.cluster_events(transformed[channels], channels, umap=False)
    clusters = adata.obs["cluster"].to_numpy()

    if len(events) > MAX_EMBED_POINTS:
        rng = np.random.default_rng(0)
        keep = np.sort(rng.choice(len(events), MAX_EMBED_POINTS, replace=False))
    else:
        keep = np.arange(len(events))
    embed = cluster.embed_events(
        transformed.iloc[keep].reset_index(drop=True), channels
    )
    cluster_ms = round((time.time() - t_cluster) * 1000)

    # what defines each cluster: per-channel median arcsinh intensity,
    # so the report answers "which markers are these" not just "how many"
    labeled = transformed[channels].assign(cluster=clusters)
    profiles = {
        str(cl): {ch: round(float(v), 3) for ch, v in row.items()}
        for cl, row in labeled.groupby("cluster")[channels].median().iterrows()
    }

    return {
        "n_events": int(len(events)),
        "n_channels": len(channels),
        "phenotypic_channels": channels,
        "qc_only_channels": qc_only,
        "n_clusters_found": int(pd.Series(clusters).nunique()),
        "cluster_sizes": {
            str(k): int(v)
            for k, v in pd.Series(clusters).value_counts().sort_index().items()
        },
        "embedding": embed.astype(float).tolist(),
        "embedding_clusters": clusters[keep].tolist(),
        "cluster_profiles": profiles,
        "channels": channel_report,
        "acquisition_drift": drift_report,
        "timing": {"cluster_embed_ms": cluster_ms},
        "scope": (
            "descriptive QC on uploaded events; no manual gates supplied, "
            "so no agreement metrics are reported"
        ),
    }


def _load_uploaded(filename: str, body: bytes) -> pd.DataFrame:
    """Parse one uploaded file into an event table."""
    suffix = Path(filename).suffix.lower()
    if suffix == ".fcs":
        # fcsparser needs a path; uploads arrive as bytes
        with tempfile.NamedTemporaryFile(suffix=".fcs", delete=False) as tmp:
            tmp.write(body)
            tmp_path = tmp.name
        try:
            df, _ = io.load_fcs(tmp_path)
        finally:
            Path(tmp_path).unlink(missing_ok=True)
        return df
    if suffix == ".csv":
        import io as stdio

        return pd.read_csv(stdio.BytesIO(body))
    msg = f"unsupported file type: {filename} (accepted: .fcs, .csv)"
    raise ValueError(msg)


def create_app():
    """Build the FastAPI app. Imported lazily so the pipeline never pays
    the web-framework import cost."""
    from typing import List

    from fastapi import FastAPI, File, Request, UploadFile
    from fastapi.middleware.cors import CORSMiddleware
    from fastapi.responses import FileResponse, JSONResponse
    from fastapi.staticfiles import StaticFiles

    app = FastAPI(title="cytof-qc", version="0.1.0")
    app.add_middleware(
        CORSMiddleware, allow_origins=["*"], allow_methods=["*"]
    )

    try:
        _STATE_DIR.mkdir(parents=True, exist_ok=True)
    except OSError:
        pass
    _load_jobs()

    dist = Path(__file__).resolve().parent.parent / "app" / "dist"
    if dist.is_dir():
        app.mount("/assets", StaticFiles(directory=dist / "assets"), name="a")

    @app.middleware("http")
    async def log_request(request, call_next):
        start = time.time()
        response = await call_next(request)
        entry = {
            "ts": int(start),
            "path": request.url.path,
            "status": response.status_code,
            "ms": round((time.time() - start) * 1000),
        }
        n_events = getattr(request.state, "n_events", None)
        if n_events is not None:
            entry["n_events"] = n_events
        # stdout too — container logs outlive the filesystem, so the
        # record survives redeploys even without a mounted volume
        logging.getLogger("cytof_qc.requests").info(json.dumps(entry))
        try:
            with _REQUEST_LOG.open("a") as fh:
                fh.write(json.dumps(entry) + "\n")
        except OSError:
            logging.getLogger("cytof_qc.service").warning("request log write failed")
        return response

    @app.get("/health")
    def health():
        return {"ok": True, "version": "0.1.0"}

    @app.get("/api/stats")
    def stats():
        """Request-log summary — the same aggregate view the offline
        summarizer prints, served live so the instrumentation is
        inspectable, not just claimed."""
        entries = []
        try:
            with _REQUEST_LOG.open() as fh:
                for line in fh:
                    try:
                        entries.append(json.loads(line))
                    except json.JSONDecodeError:
                        continue
        except OSError:
            pass
        by_status = {}
        for e in entries:
            by_status[e.get("status")] = by_status.get(e.get("status"), 0) + 1
        lat = sorted(e["ms"] for e in entries if "ms" in e)
        nev = sorted(e["n_events"] for e in entries if "n_events" in e)

        def pct(xs, p):
            return xs[min(len(xs) - 1, int(len(xs) * p))] if xs else None

        return {
            "requests": len(entries),
            "by_status": by_status,
            "latency_ms": {"p50": pct(lat, 0.5), "p95": pct(lat, 0.95)},
            "n_events": {
                "min": nev[0] if nev else None,
                "max": nev[-1] if nev else None,
            },
        }

    @app.get("/", include_in_schema=False)
    def index():
        idx = dist / "index.html"
        if idx.exists():
            return FileResponse(idx)
        return JSONResponse({"service": "cytof-qc", "docs": "/docs"})

    @app.post("/api/analyze")
    async def analyze(
        request: Request, files: List[UploadFile] = File(...)
    ):
        # throttle before any parse work — counts attempts, not successes
        fwd = (request.headers.get("x-forwarded-for") or "").split(",")[0].strip()
        client_ip = fwd or (request.client.host if request.client else "unknown")
        if not _rate_limit_ok(client_ip):
            return JSONResponse(
                {"error": f"rate limited ({RATE_LIMIT_MAX}/hour)"},
                status_code=429,
            )
        if not files or len(files) > MAX_UPLOAD_FILES:
            return JSONResponse(
                {"error": f"upload 1–{MAX_UPLOAD_FILES} files"}, status_code=400
            )
        # parse synchronously so malformed uploads fail fast with a 400;
        # the expensive cluster+embed runs as a background job
        t_load = time.time()
        frames = []
        for f in files:
            if f.size is not None and f.size > MAX_FILE_BYTES:
                return JSONResponse(
                    {"error": f"{f.filename} exceeds the "
                     f"{MAX_FILE_BYTES // 1_000_000}MB file cap"},
                    status_code=400,
                )
            body = await f.read(MAX_FILE_BYTES + 1)
            if len(body) > MAX_FILE_BYTES:
                return JSONResponse(
                    {"error": f"{f.filename} exceeds the "
                     f"{MAX_FILE_BYTES // 1_000_000}MB file cap"},
                    status_code=400,
                )
            try:
                frames.append(_load_uploaded(f.filename or "upload", body))
            except Exception as exc:
                return JSONResponse({"error": str(exc)}, status_code=400)
        load_ms = round((time.time() - t_load) * 1000)

        events = pd.concat(frames, ignore_index=True)
        n_events = len(events)
        request.state.n_events = n_events
        if n_events > MAX_EVENTS:
            return JSONResponse(
                {"error": f"{n_events} events exceeds the {MAX_EVENTS} cap"},
                status_code=400,
            )
        if n_events < MIN_EVENTS:
            return JSONResponse(
                {"error": f"{n_events} events is too few to cluster "
                 f"(need >= {MIN_EVENTS})"},
                status_code=400,
            )

        with _JOBS_LOCK:
            _expire_hung_locked()
            running = sum(1 for j in _JOBS.values() if j["status"] == "running")
            if running >= MAX_RUNNING_JOBS:
                return JSONResponse(
                    {"error": "service busy, retry shortly"}, status_code=429
                )
            # expire finished jobs past TTL, bound the map by count,
            # and keep the durable log in sync
            _compact_jobs()
            job_id = uuid.uuid4().hex[:12]
            _JOBS[job_id] = {
                "status": "running",
                "started": time.time(),
                "n_events": n_events,
            }
            _persist_job(job_id, _JOBS[job_id])

        def _run():
            t_analyze = time.time()
            try:
                report = analyze_events(events)
                report["timing"].update({
                    "load_ms": load_ms,
                    "analyze_ms": round((time.time() - t_analyze) * 1000),
                })
                with _JOBS_LOCK:
                    job = _JOBS.get(job_id)
                    if job is not None:
                        job.update(status="done", report=report)
                        _persist_job(job_id, job)
            except Exception as exc:
                with _JOBS_LOCK:
                    job = _JOBS.get(job_id)
                    if job is not None:
                        job.update(status="error", error=str(exc))
                        _persist_job(job_id, job)

        threading.Thread(target=_run, daemon=True).start()
        return {"job_id": job_id, "status": "running", "n_events": n_events}

    @app.get("/api/jobs/{job_id}")
    def job_status(job_id: str):
        with _JOBS_LOCK:
            _expire_hung_locked()
            job = dict(_JOBS.get(job_id) or {})
        if not job:
            return JSONResponse({"error": "unknown job"}, status_code=404)
        if job["status"] == "done":
            return {**job["report"], "status": "done"}
        if job["status"] == "error":
            return JSONResponse(
                {"status": "error", "error": job["error"]}, status_code=200
            )
        return {
            "status": "running",
            "n_events": job["n_events"],
            "elapsed_s": round(time.time() - job["started"], 1),
        }

    # SPA static root LAST so /api/* and /health routes match first;
    # serves index.html, /assets/*, and top-level files like
    # example_report.json.
    if dist.is_dir():
        app.mount("/", StaticFiles(directory=dist, html=True), name="spa")

    return app


app = create_app()
