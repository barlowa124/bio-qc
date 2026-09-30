"""Unit tests for cultivated meat multi-omic pipeline core modules."""
import json
import pickle
import sys
from pathlib import Path

import numpy as np
import pytest

# Add project root to path
PROJ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJ / "api"))
sys.path.insert(0, str(PROJ / "notebooks"))

API_DIR = PROJ / "api"
OUT = PROJ / "p2_state_map/output"


def test_api_artifacts_exist():
    assert (API_DIR / "model.pkl").exists()
    assert (API_DIR / "scaler.pkl").exists()
    assert (API_DIR / "model_metadata.json").exists()
    assert (API_DIR / "app.py").exists()


def test_model_metadata_structure():
    meta = json.loads((API_DIR / "model_metadata.json").read_text())
    assert "genes" in meta
    assert "classes" in meta
    assert len(meta["genes"]) == 30
    assert len(meta["classes"]) == 3


def test_model_load_and_predict():
    model = pickle.load(open(API_DIR / "model.pkl", "rb"))
    scaler = pickle.load(open(API_DIR / "scaler.pkl", "rb"))
    meta = json.load(open(API_DIR / "model_metadata.json"))
    n_genes = len(meta["genes"])
    X = np.random.lognormal(0, 1, (5, n_genes)).astype(np.float32)
    Xs = scaler.transform(X)
    preds = model.predict(Xs)
    probs = model.predict_proba(Xs)
    assert preds.shape == (5,)
    assert probs.shape == (5, 3)
    assert np.allclose(probs.sum(axis=1), 1.0)


def test_shap_results_exist():
    assert (OUT / "shap_dnn_results.json").exists()
    data = json.loads((OUT / "shap_dnn_results.json").read_text())
    assert "shap" in data
    assert "dnn" in data
    assert "cv_mean" in data["dnn"]
    assert data["dnn"]["cv_mean"] > 0.8


def test_literature_results_exist():
    assert (OUT / "literature_drug_panel_noise.json").exists()
    data = json.loads((OUT / "literature_drug_panel_noise.json").read_text())
    assert "literature_benchmark" in data or "literature_overlap" in data
    assert "drug_predictions" in data
    assert "minimal_panel" in data
    assert "qpcr_noise_simulation" in data or "noise_simulation" in data


def test_tf_ppi_results_exist():
    assert (OUT / "tf_ppi_results.json").exists()
    data = json.loads((OUT / "tf_ppi_results.json").read_text())
    assert "tf_enrichment" in data
    assert "ppi_network" in data


def test_grant_proposal_exists():
    assert (OUT / "grant_proposal_draft.md").exists()
    text = (OUT / "grant_proposal_draft.md").read_text()
    assert "NIH" in text or "R21" in text
    assert "$" in text


def test_vae_bayesian_results_exist():
    assert (OUT / "vae_bayesian_results.json").exists()
    data = json.loads((OUT / "vae_bayesian_results.json").read_text())
    assert "vae" in data
    assert "bayesian_gmm" in data or "bayesian" in data


def test_wgcna_results_exist():
    assert (OUT / "wgcna_de_batch.json").exists()
    data = json.loads((OUT / "wgcna_de_batch.json").read_text())
    assert "modules" in data or "wgcna" in data


def test_bootstrap_results_exist():
    assert (OUT / "bootstrap_ml_timeseries.json").exists()
    data = json.loads((OUT / "bootstrap_ml_timeseries.json").read_text())
    assert "bootstrap" in data
    bootstrap = data["bootstrap"]
    acc = bootstrap.get("cv_mean", bootstrap.get("mean_accuracy", 0))
    assert acc > 0.90


def test_cross_species_results_exist():
    assert (OUT / "cross_species_comparison.json").exists()


def test_cellcom_results_exist():
    assert (OUT / "cellcom_benchmark_pipeline.json").exists()


def test_state_map_exists():
    assert (OUT / "state_map_results.json").exists()
    data = json.loads((OUT / "state_map_results.json").read_text())
    assert "n_samples" in data
    assert data["n_samples"] > 0


def test_qc_panel_exists():
    assert (PROJ / "p3_qc_panel/output/qc_panel_239.json").exists()


def test_readme_exists():
    assert (PROJ / "README.md").exists()
    text = (PROJ / "README.md").read_text()
    assert "cultivated meat" in text.lower()


def test_license_exists():
    assert (PROJ / "LICENSE").exists()


def test_requirements_exists():
    assert (PROJ / "requirements.txt").exists()


def test_dockerfile_exists():
    assert (PROJ / "Dockerfile").exists()


def test_zenodo_exists():
    assert (PROJ / ".zenodo.json").exists()


@pytest.mark.skipif(not (API_DIR / "app.py").exists(), reason="API not built")
def test_flask_app_import():
    import importlib.util
    spec = importlib.util.spec_from_file_location("app", API_DIR / "app.py")
    app_module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(app_module)
    assert hasattr(app_module, "app")
