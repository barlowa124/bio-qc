"""HTTP surface for the QC pipeline — upload events, get a report.

``analyze_events`` is the product core: any caller that can produce an
event table (uploaded FCS files, CSV exports, FlowJo workspaces) gets
channel QC, a drift gate, Leiden clusters, and a UMAP embedding back.
Uploaded data carries no manual gates, so the report is descriptive —
cluster structure plus QC flags — not agreement scoring.
"""

import json
import logging
import tempfile
import time
from pathlib import Path

import numpy as np
import pandas as pd

from . import io, qc, transform

MAX_UPLOAD_FILES = 200
MAX_EVENTS = 500_000
MAX_EMBED_POINTS = 20_000

_REQUEST_LOG = Path(__file__).resolve().parent.parent / "service_requests.jsonl"


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

    transformed = transform.arcsinh(numeric)
    channel_report = qc.channel_qc(
        events.assign(**{c: numeric[c] for c in numeric.columns}),
        all_channels,
    )
    drift_report = qc.acquisition_qc(transformed, all_channels)

    from . import cluster

    adata = cluster.cluster_events(transformed[channels], channels)
    clusters = adata.obs["cluster"].to_numpy()
    embed = adata.obsm["X_umap"]

    # the browser renders a deterministic subsample; full clusters are
    # reported in the counts table either way
    if len(embed) > MAX_EMBED_POINTS:
        rng = np.random.default_rng(0)
        keep = np.sort(rng.choice(len(embed), MAX_EMBED_POINTS, replace=False))
    else:
        keep = np.arange(len(embed))

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
        "embedding": embed[keep].astype(float).tolist(),
        "embedding_clusters": clusters[keep].tolist(),
        "channels": channel_report,
        "acquisition_drift": drift_report,
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

    from fastapi import FastAPI, File, UploadFile
    from fastapi.middleware.cors import CORSMiddleware
    from fastapi.responses import FileResponse, JSONResponse
    from fastapi.staticfiles import StaticFiles

    app = FastAPI(title="cytof-qc", version="0.1.0")
    app.add_middleware(
        CORSMiddleware, allow_origins=["*"], allow_methods=["*"]
    )

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
        try:
            with _REQUEST_LOG.open("a") as fh:
                fh.write(json.dumps(entry) + "\n")
        except OSError:
            logging.getLogger("cytof_qc.service").warning("request log write failed")
        return response

    @app.get("/health")
    def health():
        return {"ok": True, "version": "0.1.0"}

    @app.get("/", include_in_schema=False)
    def index():
        idx = dist / "index.html"
        if idx.exists():
            return FileResponse(idx)
        return JSONResponse({"service": "cytof-qc", "docs": "/docs"})

    @app.post("/api/analyze")
    async def analyze(files: List[UploadFile] = File(...)):
        if not files or len(files) > MAX_UPLOAD_FILES:
            return JSONResponse(
                {"error": f"upload 1–{MAX_UPLOAD_FILES} files"}, status_code=400
            )
        frames = []
        for f in files:
            body = await f.read()
            try:
                frames.append(_load_uploaded(f.filename or "upload", body))
            except Exception as exc:
                return JSONResponse({"error": str(exc)}, status_code=400)
        try:
            return analyze_events(pd.concat(frames, ignore_index=True))
        except ValueError as exc:
            return JSONResponse({"error": str(exc)}, status_code=400)

    return app


app = create_app()
