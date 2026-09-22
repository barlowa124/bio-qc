"""ci_fixtures.py -- Generate minimal fixture files for CI testing.

This script creates dummy versions of all artifacts that the test suite
expects, without requiring the large raw data files used by the full
notebooks.  Run this in GitHub Actions before pytest.
"""
import json
import pickle
from pathlib import Path

import numpy as np

np.random.seed(42)
PROJ = Path(__file__).resolve().parent
API_DIR = PROJ / "api"
OUT = PROJ / "p2_state_map" / "output"
QC_OUT = PROJ / "p3_qc_panel" / "output"

API_DIR.mkdir(parents=True, exist_ok=True)
OUT.mkdir(parents=True, exist_ok=True)
QC_OUT.mkdir(parents=True, exist_ok=True)

# ── 1. API artifacts ──
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

panel_genes = [
    "UXS1", "PLOD1", "MALAT1", "C1D", "KIF1B", "XKR9", "LINC00574", "UGT8",
    "PPIEL", "FCRLA", "IFITM3", "PLXNA1", "UPK1B", "C11orf63", "RAB42", "HRH4",
    "NPC2", "AEBP1", "GLG1", "C3orf72", "APOL1", "SNHG3", "EMC1", "LMNA",
    "CTSA", "TOMM7", "ZNF527", "PLEKHG4B", "MRPL32", "C1QBP"
]

X_dummy = np.random.lognormal(0, 1, (100, 30)).astype(np.float32)
y_dummy = np.random.choice(["expansion_competent", "committed", "terminal"], 100, p=[0.26, 0.36, 0.38])

scaler = StandardScaler().fit(X_dummy)
model = LogisticRegression(max_iter=2000, C=1.0, random_state=42, solver="lbfgs")
model.fit(scaler.transform(X_dummy), y_dummy)

pickle.dump(model, open(API_DIR / "model.pkl", "wb"))
pickle.dump(scaler, open(API_DIR / "scaler.pkl", "wb"))
json.dump({"genes": panel_genes, "classes": list(model.classes_)}, open(API_DIR / "model_metadata.json", "w"))

api_code = '''"""Flask API for 30-gene QC panel state prediction."""
import json, pickle, numpy as np
from pathlib import Path
from flask import Flask, request, jsonify

app = Flask(__name__)
API_DIR = Path(__file__).parent
model = pickle.load(open(API_DIR / "model.pkl", "rb"))
scaler = pickle.load(open(API_DIR / "scaler.pkl", "rb"))
meta = json.load(open(API_DIR / "model_metadata.json"))

@app.route("/")
def index():
    return jsonify({"service": "Cultivated Meat 30-Gene QC Panel", "genes": meta["genes"], "states": meta["classes"], "endpoints": {"/predict": "POST gene expression values", "/health": "GET health check"}})

@app.route("/health")
def health():
    return jsonify({"status": "healthy", "model_accuracy": 0.967})

@app.route("/predict", methods=["POST"])
def predict():
    data = request.get_json()
    if not data or "expression" not in data:
        return jsonify({"error": "Missing 'expression' field."}), 400
    expr = np.array(data["expression"], dtype=np.float32).reshape(1, -1)
    if expr.shape[1] != len(meta["genes"]):
        return jsonify({"error": f"Expected {len(meta['genes'])} values, got {expr.shape[1]}"}), 400
    expr_scaled = scaler.transform(expr)
    probs = model.predict_proba(expr_scaled)[0]
    pred = model.classes_[probs.argmax()]
    return jsonify({"prediction": pred, "probabilities": {model.classes_[i]: float(p) for i, p in enumerate(probs)}, "confidence": float(probs.max())})

if __name__ == "__main__":
    app.run(host="0.0.0.0", port=5000, debug=False)
'''
(API_DIR / "app.py").write_text(api_code)

# ── 2. Output JSON fixtures ──
# shap_dnn_results.json
json.dump({
    "shap": {"top_genes": ["LMNA", "PLOD1", "C1QBP"]},
    "dnn": {"cv_mean": 0.967, "cv_std": 0.012}
}, open(OUT / "shap_dnn_results.json", "w"))

# literature_drug_panel_noise.json
json.dump({
    "literature_overlap": {"count": 15},
    "drug_predictions": [{"compound": "Rapamycin", "score": 0.85, "mechanism": "mTOR inhibition"}],
    "minimal_panel": {"n_genes": 25},
    "qpcr_noise_simulation": [{"cv_level": 0.05, "mean_accuracy": 0.95}]
}, open(OUT / "literature_drug_panel_noise.json", "w"))

# tf_ppi_results.json
json.dump({
    "tf_enrichment": {"top_tfs": ["SP1", "MYOD1"]},
    "ppi_network": {"source": "STRING", "genes_queried": panel_genes, "interactions": [{"source": "C1QBP", "target": "LMNA"}], "hub_genes": ["C1QBP", "LMNA"], "n_edges": 74, "n_nodes": 30}
}, open(OUT / "tf_ppi_results.json", "w"))

# grant_proposal_draft.md
(OUT / "grant_proposal_draft.md").write_text("# Grant Proposal Draft\n\nNIH R21 application for cultivated meat QC panel.\n\nBudget: $275,000 over 2 years.")

# vae_bayesian_results.json
json.dump({
    "vae": {"latent_dim": 10, "recon_loss": 0.45},
    "bayesian_gmm": {"n_components": 3, "bic": -1200}
}, open(OUT / "vae_bayesian_results.json", "w"))

# wgcna_de_batch.json
json.dump({
    "wgcna": {"modules": {"blue": ["LMNA", "PLOD1"], "turquoise": ["C1QBP", "CTSA"]}, "n_modules": 2}
}, open(OUT / "wgcna_de_batch.json", "w"))

# bootstrap_ml_timeseries.json
json.dump({
    "bootstrap": {"cv_mean": 0.967, "cv_std": 0.01, "n_iterations": 100}
}, open(OUT / "bootstrap_ml_timeseries.json", "w"))

# cross_species_comparison.json
json.dump({
    "bovine_vs_human": {"accuracy": 0.92, "n_genes": 30},
    "porcine_vs_human": {"accuracy": 0.89, "n_genes": 30}
}, open(OUT / "cross_species_comparison.json", "w"))

# cellcom_benchmark_pipeline.json
json.dump({
    "cellcom": {"method": "CellChat", "n_interactions": 45, "top_pathways": ["NOTCH", "WNT"]}
}, open(OUT / "cellcom_benchmark_pipeline.json", "w"))

# state_map_results.json
json.dump({
    "n_samples": 239,
    "states": {"expansion_competent": 62, "committed": 86, "terminal": 91},
    "accuracy": 0.967
}, open(OUT / "state_map_results.json", "w"))

# qc_panel_239.json
json.dump({
    "panel_genes": [{"gene": g, "importance": np.random.random()} for g in panel_genes],
    "n_genes": 30,
    "method": "variance_ranking"
}, open(QC_OUT / "qc_panel_239.json", "w"))

# .zenodo.json
json.dump({
    "title": "Multi-omic state-map pipeline (methods demonstration)",
    "creators": [{"name": "barlowa124"}],
    "upload_type": "software"
}, open(PROJ / ".zenodo.json", "w"))

print("CI fixtures generated successfully.")
print(f"  API artifacts: {API_DIR}")
print(f"  Output fixtures: {OUT}")
print(f"  QC panel fixture: {QC_OUT}")
