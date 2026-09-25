"""Streamlit dashboard for the 30-gene QC panel state prediction."""
import json
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

# Setup paths
PROJ = Path(__file__).resolve().parents[1]
API_DIR = PROJ / "api"
OUT = PROJ / "p2_state_map/output"

st.set_page_config(
    page_title="Cultivated Meat QC Panel",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.title("🧬 Cultivated Meat 30-Gene QC Panel Dashboard")
st.markdown("Manufacturing-readiness state prediction for cultivated muscle tissue")

# ── Load artifacts ──
@st.cache_resource
def load_model():
    model = pickle.load(open(API_DIR / "model.pkl", "rb"))
    scaler = pickle.load(open(API_DIR / "scaler.pkl", "rb"))
    meta = json.load(open(API_DIR / "model_metadata.json"))
    return model, scaler, meta

@st.cache_data
def load_results():
    results = {}
    for f in ["shap_dnn_results.json", "literature_drug_panel_noise.json",
              "vae_bayesian_results.json", "bootstrap_ml_timeseries.json",
              "state_map_results.json", "cross_species_comparison.json"]:
        p = OUT / f
        if p.exists():
            results[f.replace(".json", "")] = json.loads(p.read_text())
    return results

try:
    model, scaler, meta = load_model()
    results = load_results()
    genes = meta["genes"]
    classes = meta["classes"]
    loaded = True
except Exception as e:
    st.error(f"Failed to load model: {e}")
    loaded = False

# ── Sidebar ──
st.sidebar.header("Navigation")
page = st.sidebar.radio("Select page", ["Prediction", "Performance", "Explainability", "Cross-Species", "Literature & Drugs", "About"])

if loaded:
    st.sidebar.markdown("---")
    st.sidebar.subheader("Model Info")
    st.sidebar.write(f"**Genes:** {len(genes)}")
    st.sidebar.write(f"**States:** {', '.join(classes)}")
    st.sidebar.write("**CV Accuracy:** 96.7%")

# ── Page: Prediction ──
if page == "Prediction" and loaded:
    st.header("🔮 State Prediction")
    st.markdown("Enter expression values (log TPM) for the 30 genes:")

    # Create input fields in 3 columns
    cols = st.columns(3)
    expr = []
    for i, g in enumerate(genes):
        with cols[i % 3]:
            val = st.number_input(g, value=0.0, step=0.1, key=f"gene_{i}")
            expr.append(val)

    if st.button("Predict", type="primary"):
        X = np.array(expr, dtype=np.float32).reshape(1, -1)
        Xs = scaler.transform(X)
        probs = model.predict_proba(Xs)[0]
        pred = model.classes_[probs.argmax()]

        col1, col2 = st.columns(2)
        with col1:
            st.subheader("Prediction")
            state_colors = {"expansion_competent": "#2ecc71", "committed": "#f39c12", "terminal": "#e74c3c"}
            color = state_colors.get(pred, "#3498db")
            st.markdown(f"<h1 style='color:{color}'>{pred.replace('_', ' ').title()}</h1>", unsafe_allow_html=True)
            st.metric("Confidence", f"{probs.max():.1%}")
        with col2:
            st.subheader("Probabilities")
            prob_df = pd.DataFrame({
                "State": [c.replace("_", " ").title() for c in model.classes_],
                "Probability": probs
            })
            fig = px.bar(prob_df, x="State", y="Probability", color="State",
                         color_discrete_map={c.replace("_", " ").title(): state_colors.get(c, "#3498db") for c in model.classes_})
            fig.update_layout(showlegend=False, yaxis_range=[0, 1])
            st.plotly_chart(fig, use_container_width=True)

    # Batch upload
    st.markdown("---")
    st.subheader("Batch Upload")
    st.markdown("Upload a CSV with columns matching the 30 gene names. Each row = one sample.")
    uploaded = st.file_uploader("Upload CSV", type="csv")
    if uploaded is not None:
        df = pd.read_csv(uploaded)
        common = [g for g in genes if g in df.columns]
        if len(common) < len(genes):
            st.warning(f"Only {len(common)}/{len(genes)} genes matched. Missing: {set(genes) - set(common)}")
            st.error("Prediction requires all 30 panel genes; fix the CSV columns and re-upload.")
        else:
            X = df[common].values.astype(np.float32)
            Xs = scaler.transform(X)
            probs = model.predict_proba(Xs)
            preds = model.classes_[probs.argmax(axis=1)]
            conf = probs.max(axis=1)
            df["prediction"] = preds
            df["confidence"] = conf
            st.dataframe(df)
            st.download_button("Download results", df.to_csv(index=False), "predictions.csv", "text/csv")

# ── Page: Performance ──
elif page == "Performance" and loaded:
    st.header("📊 Model Performance")
    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("CV Accuracy", "96.7%", "± 1.0%")
    with col2:
        st.metric("Bootstrap", "96.3%", "± 0.8%")
    with col3:
        st.metric("DNN CV", "96.2%", "± 1.5%")

    # ML comparison
    if "bootstrap_ml_timeseries" in results:
        ml = results["bootstrap_ml_timeseries"].get("ml_comparison", {})
        if ml:
            names = list(ml.keys())
            scores = [v.get("mean", 0) for v in ml.values()]
            stds = [v.get("std", 0) for v in ml.values()]
            fig = go.Figure(data=[
                go.Bar(name="Mean CV", x=names, y=scores, error_y=dict(type="data", array=stds))
            ])
            fig.update_layout(title="ML Model Comparison", yaxis_title="Accuracy", yaxis_range=[0.85, 1.0])
            st.plotly_chart(fig, use_container_width=True)

    # Noise robustness
    if "literature_drug_panel_noise" in results:
        noise = results["literature_drug_panel_noise"].get("noise_simulation", {})
        if noise:
            cvs = list(noise.keys())
            accs = [noise[cv]["accuracy"] for cv in cvs]
            fig = px.line(x=cvs, y=accs, markers=True, labels={"x": "CV (%)", "y": "Accuracy"}, title="qPCR Noise Robustness")
            fig.update_layout(yaxis_range=[0.8, 1.0])
            st.plotly_chart(fig, use_container_width=True)

# ── Page: Explainability ──
elif page == "Explainability" and loaded:
    st.header("🔍 SHAP Explainability")
    if "shap_dnn_results" in results:
        shap = results["shap_dnn_results"]["shap"]
        top = shap.get("top_genes", [])
        if top:
            df = pd.DataFrame(top)
            fig = px.bar(df, x="importance", y="gene", orientation="h", title="Gene Importance (SHAP)")
            st.plotly_chart(fig, use_container_width=True)

        per_state = shap.get("per_state_expression", {})
        if per_state:
            state = st.selectbox("Select state", list(per_state.keys()), format_func=lambda x: x.replace("_", " ").title())
            expr = per_state[state]
            df = pd.DataFrame({"gene": list(expr.keys()), "mean_expression": list(expr.values())})
            df = df.sort_values("mean_expression", ascending=False).head(15)
            fig = px.bar(df, x="mean_expression", y="gene", orientation="h", title=f"Top Genes in {state.replace('_', ' ').title()}")
            st.plotly_chart(fig, use_container_width=True)

# ── Page: Cross-Species ──
elif page == "Cross-Species" and loaded:
    st.header("🐄 Cross-Species Validation")
    if "cross_species_comparison" in results:
        cs = results["cross_species_comparison"]
        st.json(cs, expanded=False)

    st.markdown("""
    **Bovine** (GSE173199, 38 samples): 3-state conserved
    **Porcine**: Muscle progenitor validation
    **Human** (GSE240556, 17,541 nuclei): 21/30 panel genes detected
    """)

# ── Page: Literature & Drugs ──
elif page == "Literature & Drugs" and loaded:
    st.header("📚 Literature Benchmark & Drug Predictions")

    if "literature_drug_panel_noise" in results:
        data = results["literature_drug_panel_noise"]

        lit = data.get("literature_overlap", {})
        if lit:
            st.subheader("Literature Overlap")
            sets = list(lit.keys())
            overlaps = [lit[s].get("overlap_percent", 0) for s in sets]
            fig = px.bar(x=sets, y=overlaps, labels={"x": "Gene Set", "y": "Overlap (%)"}, title="Panel vs Canonical Gene Sets")
            st.plotly_chart(fig, use_container_width=True)

        drugs = data.get("drug_predictions", {})
        if drugs:
            st.subheader("Drug/Compound Predictions")
            compounds = list(drugs.keys())
            scores = [drugs[c].get("score", 0) for c in compounds]
            mechanisms = [drugs[c].get("mechanism", "unknown") for c in compounds]
            df = pd.DataFrame({"compound": compounds, "score": scores, "mechanism": mechanisms})
            df = df.sort_values("score", ascending=False)
            fig = px.bar(df, x="score", y="compound", color="mechanism", orientation="h", title="LINCS/Connectivity Map Predictions")
            st.plotly_chart(fig, use_container_width=True)

        minimal = data.get("minimal_panel", {})
        if minimal:
            st.subheader("Minimal Panel")
            st.write("**10 genes achieve 96.6% accuracy** — matching the full 30-gene panel.")
            st.write("Selected genes:", minimal.get("selected_genes", []))

# ── Page: About ──
elif page == "About":
    st.header("ℹ️ About")
    st.markdown("""
    This dashboard accompanies the manuscript:
    **"A 30-gene qPCR panel predicts manufacturing-readiness states in cultivated muscle tissue"**

    - **Repository:** [github.com/barlowa124/cultivated-meat-multiomic](https://github.com/barlowa124/cultivated-meat-multiomic)
    - **License:** MIT
    - **Panel cost:** $50–100 / batch
    - **Turnaround:** 4–6 hours

    ### Core Technologies
    - Logistic regression + DNN classifier
    - SHAP explainability
    - Bayesian GMM probabilistic modeling
    - Cross-species validation (bovine, porcine, human)
    - LINCS L1000 drug prediction
    """)

st.sidebar.markdown("---")
st.sidebar.caption("© 2026 Rao Lab — MIT License")
