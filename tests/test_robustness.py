"""Robustness battery: /predict endpoint input validation, model contract
invariants, and metadata consistency.

The API previously 500'd on non-numeric, nested, or NaN expressions;
the hardening is in api/app.py and pinned here via the Flask test client.
"""
import json
import pickle
import sys
from pathlib import Path

import numpy as np
import pytest

PROJ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJ / "api"))

flask = pytest.importorskip("flask")
from app import app  # noqa: E402


@pytest.fixture
def client():
    return app.test_client()


@pytest.fixture
def meta():
    return json.loads((PROJ / "api" / "model_metadata.json").read_text())


class TestPredictValidation:
    """Every malformed payload must be a clean JSON 400, never a 500."""

    N = 30

    def test_missing_field(self, client):
        assert client.post("/predict", json={}).status_code == 400

    def test_non_list_expression(self, client):
        for bad in ["not a list", 42, {"bust": 88}, "a" * self.N]:
            r = client.post("/predict", json={"expression": bad})
            assert r.status_code == 400, bad
            assert "error" in r.get_json()

    def test_wrong_length(self, client):
        r = client.post("/predict", json={"expression": [1.0] * (self.N - 1)})
        assert r.status_code == 400 and "Expected 30" in r.get_json()["error"]

    def test_non_numeric_elements(self, client):
        r = client.post("/predict", json={"expression": ["x"] * self.N})
        assert r.status_code == 400
        r = client.post("/predict", json={"expression": [["x"] * self.N]})
        assert r.status_code == 400

    def test_non_finite_rejected(self, client):
        for fill in (float("nan"), float("inf"), float("-inf")):
            r = client.post("/predict", json={"expression": [fill] * self.N})
            assert r.status_code == 400

    def test_none_payload(self, client):
        r = client.post("/predict", data="not json",
                        content_type="application/json")
        assert r.status_code in (400, 415)

    def test_valid_request_contract(self, client, meta):
        r = client.post("/predict",
                        json={"expression": list(np.random.default_rng(0)
                                             .lognormal(0, 1, 30))})
        assert r.status_code == 200
        body = r.get_json()
        assert body["prediction"] in meta["classes"]
        assert set(body["probabilities"]) == set(meta["classes"])
        assert abs(sum(body["probabilities"].values()) - 1.0) < 1e-6
        assert 0.0 <= body["confidence"] <= 1.0
        assert body["confidence"] == pytest.approx(
            max(body["probabilities"].values()))


class TestModelContract:
    """The served model, scaler, and metadata must agree — a model whose
    output classes drift from the metadata silently mislabels states."""

    def test_artifacts_load(self, meta):
        model = pickle.load(open(PROJ / "api" / "model.pkl", "rb"))
        scaler = pickle.load(open(PROJ / "api" / "scaler.pkl", "rb"))
        assert model.n_features_in_ == len(meta["genes"]) == 30
        assert scaler.n_features_in_ == len(meta["genes"])
        assert set(model.classes_) == set(meta["classes"])

    def test_probabilities_normalize(self, meta):
        model = pickle.load(open(PROJ / "api" / "model.pkl", "rb"))
        scaler = pickle.load(open(PROJ / "api" / "scaler.pkl", "rb"))
        X = np.random.default_rng(1).lognormal(0, 1, (8, 30)).astype(np.float32)
        probs = model.predict_proba(scaler.transform(X))
        assert np.allclose(probs.sum(axis=1), 1.0)
        assert np.isfinite(probs).all()

    def test_metadata_gene_panel_integrity(self, meta):
        assert len(meta["genes"]) == len(set(meta["genes"])) == 30
        assert len(meta["classes"]) == 3

    def test_extreme_magnitude_inputs_dont_crash(self, client, meta):
        # lognormal outliers / zeros — biologically odd but must not 500
        for expr in ([0.0] * 30, [1e6] * 30, [1e-9] * 30):
            r = client.post("/predict", json={"expression": expr})
            assert r.status_code == 200


class TestHealthEndpoint:
    def test_health_shape(self, client):
        body = client.get("/health").get_json()
        assert body["status"] == "healthy"

    def test_index_lists_genes(self, client, meta):
        body = client.get("/").get_json()
        assert body["genes"] == meta["genes"]
        assert body["states"] == meta["classes"]
