"""Flask API for 30-gene QC panel state prediction."""
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
    raw = data["expression"]
    if not isinstance(raw, list):
        return jsonify({"error": "'expression' must be a list of numbers."}), 400
    try:
        expr = np.array(raw, dtype=np.float32).reshape(1, -1)
    except (ValueError, TypeError):
        return jsonify({"error": "'expression' must contain only numeric values."}), 400
    if expr.shape[1] != len(meta["genes"]):
        return jsonify({"error": f"Expected {len(meta['genes'])} values, got {expr.shape[1]}"}), 400
    if not np.isfinite(expr).all():
        return jsonify({"error": "'expression' must contain only finite values."}), 400
    expr_scaled = scaler.transform(expr)
    probs = model.predict_proba(expr_scaled)[0]
    pred = model.classes_[probs.argmax()]
    return jsonify({"prediction": pred, "probabilities": {model.classes_[i]: float(p) for i, p in enumerate(probs)}, "confidence": float(probs.max())})

if __name__ == "__main__":
    app.run(host="127.0.0.1", port=5000, debug=False)
