"""Sensitivity and ablation study: normalization, feature selectors, classifiers."""
import json
import warnings
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
from sklearn.feature_selection import (
    RFE,
    SelectFromModel,
    SelectKBest,
    f_classif,
    mutual_info_classif,
)
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.naive_bayes import GaussianNB
from sklearn.neighbors import KNeighborsClassifier
from sklearn.preprocessing import (
    MinMaxScaler,
    QuantileTransformer,
    RobustScaler,
    StandardScaler,
)
from sklearn.svm import SVC

warnings.filterwarnings("ignore")
PROJ = Path(__file__).resolve().parents[1]
OUT = PROJ / "p2_state_map/output"
OUT.mkdir(parents=True, exist_ok=True)

print("=" * 60)
print("SENSITIVITY & ABLATION STUDY")
print("=" * 60)

# ── Load data ──
tpm_df = pd.read_csv(r"C:\Users\asdf\CascadeProjects\rao_lab_ml\nn_pipeline\combined_data\metaflux_scfea_pseudobulk\scfea_pseudobulk_tpm.csv", index_col=0)
flux_df = pd.read_csv(r"C:\Users\asdf\CascadeProjects\rao_lab_ml\nn_pipeline\combined_data\metaflux_scfea_pseudobulk\metaflux_combined_flux_matrix_all.csv", index_col=0)
ch = sorted(set(tpm_df.columns) & set(flux_df.columns))
Xr = np.log1p(tpm_df[ch].values.T).astype(np.float32)
Xf = flux_df[ch].values.T.astype(np.float32)
gv = Xr.var(axis=0)
Xr_f = Xr[:, gv > np.percentile(gv, 25)]

# Clustering to get labels
pr = PCA(15).fit_transform(StandardScaler().fit_transform(Xr_f))
pf = PCA(15).fit_transform(StandardScaler().fit_transform(Xf))
km = KMeans(3, random_state=42, n_init=10).fit(np.hstack([pr, pf]))
ma = Xf.mean(axis=1)
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

print(f"Data: {len(ch)} samples, {len(panel_genes)} genes, states={dict(Counter(y))}")

# ── 1. Normalization comparison ──
print("\n─── 1. Normalization Methods ───")
normalizers = {
    "StandardScaler": StandardScaler(),
    "RobustScaler": RobustScaler(),
    "MinMaxScaler": MinMaxScaler(),
    "QuantileTransformer": QuantileTransformer(n_quantiles=min(1000, len(ch)), output_distribution="normal"),
    "Log1P_only": None,  # no further scaling
}
norm_results = {}
for name, scaler in normalizers.items():
    if scaler is None:
        X_scaled = X_panel
    else:
        X_scaled = scaler.fit_transform(X_panel)
    scores = cross_val_score(LogisticRegression(C=1.0, max_iter=2000, random_state=42), X_scaled, y, cv=StratifiedKFold(5, shuffle=True, random_state=42), scoring="accuracy")
    norm_results[name] = {"mean": float(scores.mean()), "std": float(scores.std())}
    print(f"  {name}: {scores.mean():.4f} ± {scores.std():.4f}")

# ── 2. Feature selection methods ──
print("\n─── 2. Feature Selection Methods ───")
X_ss = StandardScaler().fit_transform(X_panel)
selectors = {
    "f_classif_top15": SelectKBest(f_classif, k=15),
    "mutual_info_top15": SelectKBest(mutual_info_classif, k=15),
    "logistic_rfe": RFE(LogisticRegression(C=1.0, max_iter=2000, random_state=42), n_features_to_select=15),
    "rf_importance": SelectFromModel(RandomForestClassifier(n_estimators=200, random_state=42), max_features=15),
    "all_features": None,
}
sel_results = {}
for name, selector in selectors.items():
    if selector is None:
        X_sel = X_ss
        n_feat = X_ss.shape[1]
    else:
        X_sel = selector.fit_transform(X_ss, y)
        n_feat = X_sel.shape[1]
    scores = cross_val_score(LogisticRegression(C=1.0, max_iter=2000, random_state=42), X_sel, y, cv=StratifiedKFold(5, shuffle=True, random_state=42), scoring="accuracy")
    sel_results[name] = {"mean": float(scores.mean()), "std": float(scores.std()), "n_features": n_feat}
    print(f"  {name} ({n_feat} feats): {scores.mean():.4f} ± {scores.std():.4f}")

# ── 3. Classifier comparison (more extensive) ──
print("\n─── 3. Classifier Comparison ───")
classifiers = {
    "LogisticRegression": LogisticRegression(C=1.0, max_iter=2000, random_state=42),
    "LogisticRegression_L2_001": LogisticRegression(C=0.01, penalty="l2", max_iter=2000, random_state=42),
    "LogisticRegression_L1": LogisticRegression(C=1.0, penalty="l1", solver="saga", max_iter=5000, random_state=42),
    "SVM_RBF": SVC(kernel="rbf", C=1.0, probability=True, random_state=42),
    "SVM_Linear": SVC(kernel="linear", C=1.0, probability=True, random_state=42),
    "RandomForest_200": RandomForestClassifier(n_estimators=200, max_depth=10, random_state=42),
    "RandomForest_500": RandomForestClassifier(n_estimators=500, max_depth=None, random_state=42),
    "GradientBoosting": GradientBoostingClassifier(n_estimators=200, max_depth=5, random_state=42),
    "KNN_5": KNeighborsClassifier(n_neighbors=5),
    "KNN_10": KNeighborsClassifier(n_neighbors=10),
    "NaiveBayes": GaussianNB(),
}
clf_results = {}
for name, clf in classifiers.items():
    scores = cross_val_score(clf, X_ss, y, cv=StratifiedKFold(5, shuffle=True, random_state=42), scoring="accuracy")
    clf_results[name] = {"mean": float(scores.mean()), "std": float(scores.std())}
    print(f"  {name}: {scores.mean():.4f} ± {scores.std():.4f}")

# ── 4. Ablation: remove each gene ──
print("\n─── 4. Leave-One-Gene-Out Ablation ───")
base_scores = cross_val_score(LogisticRegression(C=1.0, max_iter=2000, random_state=42), X_ss, y, cv=StratifiedKFold(5, shuffle=True, random_state=42), scoring="accuracy")
base_mean = base_scores.mean()
ablation = {}
for i, g in enumerate(panel_genes):
    X_loo = np.delete(X_ss, i, axis=1)
    scores = cross_val_score(LogisticRegression(C=1.0, max_iter=2000, random_state=42), X_loo, y, cv=StratifiedKFold(5, shuffle=True, random_state=42), scoring="accuracy")
    delta = base_mean - scores.mean()
    ablation[g] = {"without_accuracy": float(scores.mean()), "delta": float(delta)}
    print(f"  Without {g}: {scores.mean():.4f} (Δ={delta:.4f})")

# Sort by importance (largest delta = most important)
ablation_sorted = dict(sorted(ablation.items(), key=lambda x: x[1]["delta"], reverse=True))

# ── 5. Hyperparameter grid for LogisticRegression ──
print("\n─── 5. Logistic Regression Hyperparameter Sensitivity ───")
C_values = [0.001, 0.01, 0.1, 1.0, 10.0, 100.0]
hyperparam_results = {}
for C in C_values:
    scores = cross_val_score(LogisticRegression(C=C, max_iter=2000, random_state=42), X_ss, y, cv=StratifiedKFold(5, shuffle=True, random_state=42), scoring="accuracy")
    hyperparam_results[str(C)] = {"mean": float(scores.mean()), "std": float(scores.std())}
    print(f"  C={C}: {scores.mean():.4f} ± {scores.std():.4f}")

# ── Save ──
results = {
    "normalization": norm_results,
    "feature_selection": sel_results,
    "classifiers": clf_results,
    "ablation": ablation_sorted,
    "hyperparameters": hyperparam_results,
    "base_accuracy": float(base_mean),
}
json.dump(results, open(OUT / "sensitivity_ablation.json", "w"), indent=2)
print("\nSaved to sensitivity_ablation.json")
print("DONE")
