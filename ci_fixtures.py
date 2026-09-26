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
    app.run(host="127.0.0.1", port=5000, debug=False)
'''
(API_DIR / "app.py").write_text(api_code)

# ── 2. Output JSON fixtures ──
# shap_dnn_results.json (values per audited fig2 CV benchmarks)
json.dump({
    "data_source": "ci_fixture_synthetic",
    "shap": {"top_genes": [
        {"gene": "LMNA", "importance": 0.09115},
        {"gene": "PLOD1", "importance": 0.05773},
        {"gene": "C1QBP", "importance": 0.0727}
    ]},
    "dnn": {"cv_mean": 0.962, "cv_std": 0.015}
}, open(OUT / "shap_dnn_results.json", "w"))

# literature_drug_panel_noise.json (schema matches notebooks/literature_drug_panel.py)
json.dump({
    "data_source": "ci_fixture_synthetic",
    "literature_benchmark": [
        {"gene_set": "PanglaoDB_Satellite_Cell", "size": 20, "overlap_count": 0,
         "overlap_genes": [], "jaccard": 0.0, "pct_of_panel": 0.0, "pct_of_set": 0.0},
        {"gene_set": "Muscle_Differentiation_Canonical", "size": 20, "overlap_count": 0,
         "overlap_genes": [], "jaccard": 0.0, "pct_of_panel": 0.0, "pct_of_set": 0.0}
    ],
    "drug_predictions": [
        {"compound": "mTOR_inhibitor_Rapamycin", "target": "MTOR",
         "effect": "promotes_quiescence", "ref": "Rodgers 2014",
         "panel_gene_matches": [], "relevance_score": 0},
        {"compound": "p38_inhibitor_SB203580", "target": "MAPK14",
         "effect": "maintains_quiescence", "ref": "Bernet 2014",
         "panel_gene_matches": [], "relevance_score": 0}
    ],
    "minimal_panel": {"results": [], "knee_point_n_genes": 10},
    "qpcr_noise_simulation": [{"cv_level": 0.05, "mean_accuracy": 0.95,
                               "std_accuracy": 0.02, "min_accuracy": 0.90,
                               "accuracy_loss": 0.017}]
}, open(OUT / "literature_drug_panel_noise.json", "w"))

# tf_ppi_results.json
json.dump({
    "data_source": "ci_fixture_synthetic",
    "tf_enrichment": {"top_tfs": ["SP1", "MYOD1"]},
    "ppi_network": {"source": "STRING", "genes_queried": panel_genes, "interactions": [{"source": "C1QBP", "target": "LMNA"}], "hub_genes": ["C1QBP", "LMNA"], "n_edges": 74, "n_nodes": 30}
}, open(OUT / "tf_ppi_results.json", "w"))

# grant_proposal_draft.md
(OUT / "grant_proposal_draft.md").write_text("# Grant Proposal Draft\n\nNIH R21 application for cultivated meat QC panel.\n\nBudget: $275,000 over 2 years.")

# vae_bayesian_results.json (script constants: LATENT_DIM=8, GMM n_components=3)
json.dump({
    "data_source": "ci_fixture_synthetic",
    "vae": {"latent_dim": 8},
    "bayesian_gmm": {"n_components": 3}
}, open(OUT / "vae_bayesian_results.json", "w"))

# wgcna_de_batch.json
json.dump({
    "data_source": "ci_fixture_synthetic",
    "wgcna": {"modules": {"blue": ["LMNA", "PLOD1"], "turquoise": ["C1QBP", "CTSA"]}, "n_modules": 2}
}, open(OUT / "wgcna_de_batch.json", "w"))

# bootstrap_ml_timeseries.json (n_boot=1000 in script; 96.3% +/- 0.8% per docs)
json.dump({
    "data_source": "ci_fixture_synthetic",
    "bootstrap": {"cv_mean": 0.963, "cv_std": 0.008, "n_iterations": 1000}
}, open(OUT / "bootstrap_ml_timeseries.json", "w"))

# cross_species_comparison.json (GSE240556 is bovine snRNA-seq, not human)
json.dump({
    "data_source": "ci_fixture_synthetic",
    "bovine_vs_human": {"accuracy": 0.92, "n_genes": 30},
    "porcine_vs_human": {"accuracy": 0.88, "n_genes": 30},
    "bovine_snrna_vs_human": {"accuracy": 0.85, "n_genes": 30,
                              "n_samples": 17541, "note": "GSE240556 is Bos taurus"}
}, open(OUT / "cross_species_comparison.json", "w"))

# cellcom_benchmark_pipeline.json
json.dump({
    "data_source": "ci_fixture_synthetic",
    "cellcom": {"method": "CellChat", "n_interactions": 45, "top_pathways": ["NOTCH", "WNT"]}
}, open(OUT / "cellcom_benchmark_pipeline.json", "w"))

# state_map_results.json (cluster sizes 61/86/92 per MODEL_CARD and poster)
json.dump({
    "data_source": "ci_fixture_synthetic",
    "n_samples": 239,
    "states": {"expansion_competent": 61, "committed": 86, "terminal": 92},
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
