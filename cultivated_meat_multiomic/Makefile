"""Makefile for one-command reproducibility of the cultivated meat multi-omic pipeline."""
# Makefile for Cultivated Meat Multi-Omic Analysis Pipeline
# Run `make all` to execute the full computational pipeline

.PHONY: all setup test lint docker docs dashboard api clean

PYTHON := python
PIP := pip

# Directories
NOTEBOOKS := notebooks
API_DIR := api
DASHBOARD_DIR := dashboard
TESTS := tests
DOCS_DIR := docs

# Default target
all: setup test analysis api report

# Install dependencies
setup:
	$(PIP) install -r requirements.txt
	$(PIP) install pytest black flake8 mypy streamlit plotly

# Run unit tests
test:
	pytest $(TESTS)/ -v --tb=short

# Code quality checks
lint:
	black $(NOTEBOOKS)/ --check || black $(NOTEBOOKS)/
	flake8 $(NOTEBOOKS)/ --count --select=E9,F63,F7,F82 --show-source --statistics
	mypy $(NOTEBOOKS)/ --ignore-missing-imports || true

# Run core analyses
analysis: state_map qc_panel bootstrap shap_dnn wgcna cellcom literature tf_ppi vae_bayesian

state_map:
	$(PYTHON) $(NOTEBOOKS)/p2_state_map.py

qc_panel:
	$(PYTHON) $(NOTEBOOKS)/p3_qc_panel.py

bootstrap:
	$(PYTHON) $(NOTEBOOKS)/bootstrap_ml_timeseries.py

shap_dnn:
	$(PYTHON) $(NOTEBOOKS)/shap_dnn.py

wgcna:
	$(PYTHON) $(NOTEBOOKS)/wgcna_de_batch.py

cellcom:
	$(PYTHON) $(NOTEBOOKS)/cellcom_benchmark.py

literature:
	$(PYTHON) $(NOTEBOOKS)/literature_drug_panel.py

tf_ppi:
	$(PYTHON) $(NOTEBOOKS)/tf_ppi_analysis.py

vae_bayesian:
	$(PYTHON) $(NOTEBOOKS)/vae_bayesian.py

# Build API artifacts
api:
	$(PYTHON) $(NOTEBOOKS)/build_api.py

# Build reports
report:
	$(PYTHON) $(NOTEBOOKS)/build_report.py

# Build documentation site
docs:
	$(PYTHON) $(NOTEBOOKS)/build_docs.py

# Run dashboard
dashboard:
	streamlit run $(DASHBOARD_DIR)/app.py

# Build Docker image
docker:
	docker build -t cultivated-meat-multiomic .

# Run Docker container
docker-run:
	docker run --rm -p 5000:5000 cultivated-meat-multiomic api/app.py

# Clean generated files
clean:
	find . -type f -name "*.pyc" -delete
	find . -type d -name "__pycache__" -delete
	find . -type d -name ".pytest_cache" -exec rm -rf {} + || true
	find . -type d -name ".mypy_cache" -exec rm -rf {} + || true
