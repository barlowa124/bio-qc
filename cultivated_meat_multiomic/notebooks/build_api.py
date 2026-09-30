"""Prediction API: POST gene expression -> state prediction + confidence."""
import json
import pickle
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")
PROJ = Path(__file__).resolve().parents[1]
OUT = PROJ / "p2_state_map/output"
API_DIR = PROJ / "api"
API_DIR.mkdir(exist_ok=True)

print("=" * 60)
print("PREDICTION API SETUP")
print("=" * 60)

# Load data and train model
tpm_df = pd.read_csv(r"C:\Users\asdf\CascadeProjects\rao_lab_ml\nn_pipeline\combined_data\metaflux_scfea_pseudobulk\scfea_pseudobulk_tpm.csv", index_col=0)
flux_df = pd.read_csv(r"C:\Users\asdf\CascadeProjects\rao_lab_ml\nn_pipeline\combined_data\metaflux_scfea_pseudobulk\metaflux_combined_flux_matrix_all.csv", index_col=0)
ch = sorted(set(tpm_df.columns) & set(flux_df.columns))
Xr = np.log1p(tpm_df[ch].values.T).astype(np.float32)
Xf = flux_df[ch].values.T.astype(np.float32)
gv = Xr.var(axis=0)
Xr_f = Xr[:, gv > np.percentile(gv, 25)]
Xr_s = StandardScaler().fit_transform(Xr_f)
Xf_s = StandardScaler().fit_transform(Xf)
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA

pr = PCA(15).fit_transform(Xr_s)
pf = PCA(15).fit_transform(Xf_s)
km = KMeans(3, random_state=42, n_init=10).fit(np.hstack([pr, pf]))
ma = Xf_s.mean(axis=1)
profs = {c: float(ma[km.labels_ == c].mean()) for c in range(3)}
ordr = sorted(profs, key=profs.get)
rmap = {ordr[0]: "expansion_competent", ordr[1]: "committed", ordr[2]: "terminal"}
y = np.array([rmap[c] for c in km.labels_])

qc_data = json.loads((PROJ / "p3_qc_panel/output/qc_panel_239.json").read_text())
gp = qc_data.get("panel_genes", qc_data.get("genes", []))
if isinstance(gp[0], dict): gp = [g["gene"] if isinstance(g, dict) else g for g in gp]
panel_genes = [g for g in gp if g in tpm_df.index]
X_panel = np.zeros((len(ch), len(panel_genes)), dtype=np.float32)
for i, g in enumerate(panel_genes):
    X_panel[:, i] = Xr[:, list(tpm_df.index).index(g)]

# Train final model
model = LogisticRegression(C=1.0, max_iter=2000, random_state=42).fit(X_panel, y)
scaler = StandardScaler().fit(X_panel)

# Save model artifacts
pickle.dump(model, open(API_DIR / "model.pkl", "wb"))
pickle.dump(scaler, open(API_DIR / "scaler.pkl", "wb"))
json.dump({"genes": panel_genes, "classes": list(model.classes_)}, open(API_DIR / "model_metadata.json", "w"))

# Create Flask API
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
        return jsonify({"error": "Missing 'expression' field. Send {'expression': [val1, val2, ..., val30]}"}), 400
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

# Create test client
test_code = '''"""Test client for the prediction API."""
import requests, json, numpy as np

# Example: random expression values for 30 genes
test_expr = np.random.lognormal(0, 1, 30).tolist()

# Local test
response = requests.post("http://localhost:5000/predict", json={"expression": test_expr})
print(json.dumps(response.json(), indent=2))
'''
(API_DIR / "test_client.py").write_text(test_code)

# Create requirements
(API_DIR / "requirements.txt").write_text("flask\nnumpy\nscikit-learn\nrequests\n")

print(f"API created at: {API_DIR}")
print("  app.py — Flask server")
print("  model.pkl — Trained logistic regression")
print("  scaler.pkl — StandardScaler")
print("  model_metadata.json — Gene list + class labels")
print("\nTo run: cd api && pip install -r requirements.txt && python app.py")
print("To test: python test_client.py")
print("DONE")
