"""Deterministic synthetic organoid + reference data for the demo path.

Generates a reference "tissue" with one population per marker key in
config.markers. Each reference type expresses a private 50-gene program
(real cell types differ across hundreds of genes, not just markers) plus
its named markers at the highest rate.

The organoid culture:
  - recovers all-but-one reference population at reduced program rate
    (immature phenotype),
  - contains one aberrant cluster with its own private program that
    correlates with no reference type.

This is a test fixture with known ground truth. Never present its scores
as biological findings.
"""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd

_PROGRAM_SIZE = 50


def _sample_population(rng, centroid, features, n_cells, dropout):
    """Poisson counts on gamma-dispersed rates + Bernoulli dropout."""
    n_cells = int(n_cells)
    lam = rng.gamma(shape=2.0, scale=centroid / 2.0, size=(n_cells, len(centroid)))
    lam *= rng.binomial(1, 1 - dropout, size=lam.shape)
    for g, rate in features:
        lam[:, g] += rng.gamma(shape=2.0, scale=rate / 2.0, size=n_cells)
    return rng.poisson(lam)


def make_demo(cfg: dict) -> tuple[ad.AnnData, ad.AnnData]:
    demo = cfg["dataset"]["demo"]
    rng = np.random.default_rng(demo["seed"])
    n_genes = demo["n_genes"]
    marker_sets = cfg["dataset"]["demo"]["markers"]
    cell_types = list(marker_sets)

    # Layout: per type, 2 named marker genes then a 50-gene program;
    # then a private program for the aberrant organoid population.
    var_names = [f"GENE_{i:04d}" for i in range(n_genes)]
    marker_idx, program_idx = {}, {}
    pos = 0
    for ct, genes in marker_sets.items():
        marker_idx[ct] = []
        for g in genes:
            var_names[pos] = g
            marker_idx[ct].append(pos)
            pos += 1
        program_idx[ct] = list(range(pos, pos + _PROGRAM_SIZE))
        pos += _PROGRAM_SIZE
    aberr_program = list(range(pos, pos + _PROGRAM_SIZE))
    pos += _PROGRAM_SIZE
    if pos > n_genes:
        raise ValueError(f"demo needs >={pos} genes, got {n_genes}")

    base_rate = np.full(n_genes, 2.0)

    # Reference: balanced populations, program + markers fully on.
    ref_counts, ref_labels = [], []
    per_type = demo["n_reference_cells"] // len(cell_types)
    for ct in cell_types:
        feats = [(g, 60.0) for g in marker_idx[ct]] + [
            (g, 25.0) for g in program_idx[ct]
        ]
        ref_counts.append(
            _sample_population(rng, base_rate, feats, per_type, 0.1)
        )
        ref_labels += [ct] * per_type
    reference = ad.AnnData(
        X=np.vstack(ref_counts).astype(np.float32),
        obs=pd.DataFrame(
            {"cell_type": ref_labels},
            index=[f"ref_{i}" for i in range(sum(len(r) for r in ref_counts))],
        ),
        var=pd.DataFrame(index=var_names),
    )

    # Organoid: all-but-one type at reduced program rate (immature) plus
    # an aberrant population on its own private program.
    n_org = demo["n_organoid_cells"]
    present, missing = cell_types[:-1], cell_types[-1]
    # 75% of organoid cells are typed; the remainder is the aberrant pool
    typed_fraction = 0.75
    org_counts, org_labels = [], []
    per_type = int(n_org * typed_fraction) // len(present)
    for ct in present:
        feats = [(g, 25.0) for g in marker_idx[ct]] + [
            (g, 10.0) for g in program_idx[ct]
        ]
        org_counts.append(
            _sample_population(rng, base_rate, feats, per_type, 0.2)
        )
        org_labels += [ct] * per_type
    n_aberr = n_org - sum(len(c) for c in org_counts)
    org_counts.append(
        _sample_population(
            rng,
            base_rate,
            [(g, 30.0) for g in aberr_program],
            n_aberr,
            0.15,
        )
    )
    org_labels += ["unannotated"] * n_aberr
    organoid = ad.AnnData(
        X=np.vstack(org_counts).astype(np.float32),
        obs=pd.DataFrame(
            {
                "cell_type": org_labels,
                "true_population": org_labels,
            },
            index=[
                f"org_{i}" for i in range(sum(len(c) for c in org_counts))
            ],
        ),
        var=pd.DataFrame(index=var_names),
    )
    organoid.uns["demo_missing_type"] = missing
    organoid.uns["markers"] = {k: list(v) for k, v in marker_sets.items()}
    reference.uns["markers"] = {k: list(v) for k, v in marker_sets.items()}
    return organoid, reference
