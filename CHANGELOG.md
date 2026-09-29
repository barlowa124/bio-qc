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

### Hardening — second review pass
- Per-file size cap (`MAX_FILE_BYTES`) checked pre-read and enforced
  with a bounded read so a missing/misleading size header can't buffer
  unbounded data; file-count and total-event caps were already in.
- `MIN_EVENTS` gate: tiny uploads now 400 immediately instead of
  crashing inside neighbors/UMAP. Verified live: 5-event CSV returns
  `{"error":"5 events is too few to cluster (need >= 16)"}`.
- `JOB_TTL_S` was declared but unused — finished jobs now actually
  expire, and the registry stays bounded.
- `MAX_JOB_RUNTIME_S`: a hung worker (not a crash — exceptions were
  already caught) used to hold a concurrency slot forever; stale
  running jobs now transition to `error` on both poll and submit paths.
- Example report: `app/public/example_report.json` (9,222-event
  Levine_13dim subset, precomputed) renders via the "load an example
  report" link — reviewers can see the product with no upload.
- Fixed while shipping it: `app/public/` wasn't in the Dockerfile's
  frontend COPY list, so the file silently never reached the image;
  and only `dist/assets` was mounted, so top-level dist files 404'd.
  Both caught by testing the deployed URL, not the local build.
- `GET /api/stats` serves the request-log summary live (counts,
  status split, p50/p95 latency, event-range) — the instrumentation
  is inspectable, not just claimed. Guard paths that were previously
  untested now covered: file-size cap via bounded read, 429 at the
  concurrency limit, TTL expiry of finished jobs.
