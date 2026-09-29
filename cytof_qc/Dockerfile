# cytof-qc service: FastAPI backend + built React frontend
FROM node:20-slim AS frontend
WORKDIR /app
COPY app/package.json app/package-lock.json app/tsconfig.json app/vite.config.ts app/index.html ./
COPY app/src ./src
COPY app/public ./public
RUN npm ci && npm run build

FROM python:3.11-slim
WORKDIR /srv
COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt
COPY cytof_qc ./cytof_qc
COPY --from=frontend /app/dist ./app/dist
EXPOSE 8000
# Railway and similar platforms inject $PORT; default 8000 locally
CMD ["sh", "-c", "uvicorn cytof_qc.service:app --host 0.0.0.0 --port ${PORT:-8000}"]
