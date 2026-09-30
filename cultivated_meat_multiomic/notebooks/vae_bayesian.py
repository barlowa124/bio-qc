"""VAE embeddings + Bayesian state assignment."""
import json
import warnings
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.optim as optim
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.mixture import BayesianGaussianMixture
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")
PROJ = Path(__file__).resolve().parents[1]
OUT = PROJ / "p2_state_map/output"
OUT.mkdir(exist_ok=True)

print("=" * 60)
print("VAE EMBEDDINGS + BAYESIAN STATE ASSIGNMENT")
print("=" * 60)

# Load data
tpm_df = pd.read_csv(r"C:\Users\asdf\CascadeProjects\rao_lab_ml\nn_pipeline\combined_data\metaflux_scfea_pseudobulk\scfea_pseudobulk_tpm.csv", index_col=0)
flux_df = pd.read_csv(r"C:\Users\asdf\CascadeProjects\rao_lab_ml\nn_pipeline\combined_data\metaflux_scfea_pseudobulk\metaflux_combined_flux_matrix_all.csv", index_col=0)
ch = sorted(set(tpm_df.columns) & set(flux_df.columns))
Xr = np.log1p(tpm_df[ch].values.T).astype(np.float32)
Xf = flux_df[ch].values.T.astype(np.float32)
gv = Xr.var(axis=0)
kp = gv > np.percentile(gv, 25)
Xr_f = Xr[:, kp]
Xr_s = StandardScaler().fit_transform(Xr_f)
Xf_s = StandardScaler().fit_transform(Xf)
Xj = np.hstack([Xr_s, Xf_s]).astype(np.float32)

# ── VAE ──
class VAE(nn.Module):
    def __init__(self, in_dim, latent=8, hidden=128):
        super().__init__()
        self.enc = nn.Sequential(nn.Linear(in_dim, hidden), nn.ReLU(), nn.Linear(hidden, 64), nn.ReLU())
        self.mu = nn.Linear(64, latent)
        self.logvar = nn.Linear(64, latent)
        self.dec = nn.Sequential(nn.Linear(latent, 64), nn.ReLU(), nn.Linear(64, hidden), nn.ReLU(), nn.Linear(hidden, in_dim))

    def encode(self, x):
        h = self.enc(x)
        return self.mu(h), self.logvar(h)

    def reparam(self, mu, logvar):
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        return mu + eps * std

    def forward(self, x):
        mu, logvar = self.encode(x)
        z = self.reparam(mu, logvar)
        return self.dec(z), mu, logvar

def loss_fn(recon, x, mu, logvar):
    mse = nn.MSELoss()(recon, x)
    kl = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp())
    return mse + 0.001 * kl / x.size(0)

# Train VAE
in_dim = Xj.shape[1]
vae = VAE(in_dim, latent=8)
opt = optim.Adam(vae.parameters(), lr=0.001)
Xj_t = torch.tensor(Xj)

vae.train()
for epoch in range(200):
    opt.zero_grad()
    recon, mu, logvar = vae(Xj_t)
    loss = loss_fn(recon, Xj_t, mu, logvar)
    loss.backward()
    opt.step()
    if (epoch + 1) % 50 == 0:
        print(f"  Epoch {epoch+1}: loss={loss.item():.4f}")

# Get VAE embeddings
vae.eval()
with torch.no_grad():
    mu, _ = vae.encode(Xj_t)
    vae_emb = mu.numpy()

# Compare PCA vs VAE
pca_emb = PCA(n_components=8).fit_transform(Xj)
pca_var = PCA(n_components=8).fit(Xj).explained_variance_ratio_.sum()
vae_var = np.var(vae_emb, axis=0).sum() / np.var(Xj, axis=0).sum()
print(f"\nPCA 8D variance retained: {pca_var:.3f}")
print(f"VAE 8D variance ratio: {vae_var:.3f}")

# ── Bayesian Gaussian Mixture ──
print("\n─── Bayesian GMM State Assignment ───")
bgm = BayesianGaussianMixture(n_components=3, random_state=42, max_iter=500, weight_concentration_prior_type='dirichlet_process')
bgm.fit(vae_emb)
probs = bgm.predict_proba(vae_emb)
hard = bgm.predict(vae_emb)

# Order by metabolic activity
ma = Xf_s.mean(axis=1)
profs = {}
for c in range(3):
    m = hard == c
    profs[c] = float(ma[m].mean())
ordr = sorted(profs, key=lambda c: profs[c])
rmap = {ordr[0]: "expansion_competent", ordr[1]: "committed", ordr[2]: "terminal"}
states = np.array([rmap[c] for c in hard])

print(f"Bayesian states: {dict(Counter(states))}")

# Uncertainty quantification
max_probs = probs.max(axis=1)
confident = (max_probs > 0.8).sum()
uncertain = (max_probs <= 0.8).sum()
print(f"Confident assignments (>0.8): {confident} ({confident/len(max_probs):.1%})")
print(f"Uncertain assignments (<=0.8): {uncertain} ({uncertain/len(max_probs):.1%})")

# Per-sample probabilities
sample_probs = []
for i, s in enumerate(ch):
    p = probs[i]
    ordered = {rmap[j]: float(p[j]) for j in range(3)}
    sample_probs.append({"sample": s, "state": states[i], "max_prob": float(max_probs[i]), "probs": ordered})

# Compare with K-means
km = KMeans(n_clusters=3, random_state=42, n_init=10)
km_labels = km.fit_predict(vae_emb)
agreement = (hard == km_labels).mean()
print(f"BGM-KMeans agreement: {agreement:.3f}")

# Save
results = {
    "vae": {"latent_dim": 8, "pca_variance_retained": float(pca_var), "vae_variance_ratio": float(vae_var)},
    "bayesian_gmm": {
        "n_components": 3, "state_distribution": {s: int(v) for s, v in Counter(states).items()},
        "confident_assignments": int(confident), "uncertain_assignments": int(uncertain),
        "bgm_kmeans_agreement": float(agreement)
    },
    "sample_probabilities": sample_probs
}
json.dump(results, open(OUT / "vae_bayesian_results.json", "w"), indent=2)
print("\nSaved to vae_bayesian_results.json")
print("=" * 60)
print("DONE")
