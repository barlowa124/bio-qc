# Cultivated Meat Multi-Omic Analysis Pipeline
# Docker image for one-command reproducibility

FROM python:3.11-slim

LABEL maintainer="Rao Lab"
LABEL description="Multi-omic state map + 30-gene QC panel for cultivated meat manufacturing"

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

ENV PYTHONUNBUFFERED=1
ENV PYTHONDONTWRITEBYTECODE=1

RUN mkdir -p /app/p2_state_map/output /app/p3_qc_panel/output \
    /app/p1_protocol_meta_analysis/output /app/p4_media_formulation/output \
    /app/p5_reproducibility/output

VOLUME ["/app/data", "/app/output"]

ENTRYPOINT ["python"]
CMD ["notebooks/comprehensive_analysis.py"]
