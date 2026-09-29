# Changelog

Changes driven by measured behavior on the deployed service
(Railway, shared CPU). Timings from `service_requests.jsonl` stage
fields (`load_ms`, `cluster_embed_ms`, `analyze_ms`).

## 2026-09-28

### Added
- Web service: `POST /api/analyze` (FCS/CSV upload), React/TypeScript
  report viewer, GHCR-publishable Dockerfile, `GET /health`,
  JSONL request logging, `scripts/summarize_requests.py` for reading
  the log back.

### Changed — first instrumented iterations
- **Embed only the displayed subsample.** First real FCS upload
  (51k events, 5 Levine_13dim populations) measured `analyze_ms` at
  124s, with `load_ms` at 6ms — UMAP over all events was the cost even
  though the browser only ever receives ≤20k points. Clustering still
  runs on all events; UMAP now runs on the display subsample.
  Measured result: 124s → 90s.
- **Analysis runs as a job.** 90s is still too long to hold an HTTP
  request open, so `POST /api/analyze` parses synchronously (bad
  uploads still 400 fast) then returns `{job_id}`;
  `GET /api/jobs/{id}` polls to completion. Verified on production:
  submit returns immediately, 800-event job completes in ~15s.
