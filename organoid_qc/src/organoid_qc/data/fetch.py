"""Fetch real data from CELLxGENE Census.

Pulls a subsampled organoid slice (default: human intestine organoid cell
atlas) and a reference slice of human intestinal tissue epithelial cells,
writes both as .h5ad under data/raw/, and prints the obs paths for
config.dataset.organoid_h5ad / reference_h5ad.

Not part of the default DAG. Run once to stage real inputs:

    .venv/bin/python -m organoid_qc.data.fetch

Census access is read-only public data; nothing licensed-restricted is
downloaded beyond standard CELLxGENE terms.
"""

from __future__ import annotations

import sys

import anndata as ad
import cellxgene_census
from organoid_qc.util import load_config


def fetch(cfg: dict) -> tuple[ad.AnnData, ad.AnnData]:
    f = cfg["fetch"]
    census = cellxgene_census.open_soma(census_version=f["census_version"])
    try:
        exp = census["census_data"]["homo_sapiens"]

        # Organoid: subsample cells from the configured atlas dataset.
        obs = (
            exp.obs.read(
                value_filter=f"dataset_id == '{f['organoid_dataset_id']}'",
                column_names=["soma_joinid"],
            )
            .concat()
            .to_pandas()
        )
        org_ids = obs["soma_joinid"].sample(
            min(f["n_organoid_cells"], len(obs)),
            random_state=f["seed"],
        )
        organoid = cellxgene_census.get_anndata(
            census, "Homo sapiens", obs_coords=org_ids.tolist()
        )

        # Reference: human intestinal tissue, the panel's cell types.
        ct_filter = " or ".join(
            f"cell_type == '{ct}'" for ct in f["reference_cell_types"]
        )
        tissue_filter = " or ".join(
            f"tissue_general == '{t}'" for t in f["reference_tissues"]
        )
        obs = (
            exp.obs.read(
                value_filter=f"({ct_filter}) and ({tissue_filter})",
                column_names=["soma_joinid", "cell_type"],
            )
            .concat()
            .to_pandas()
        )
        ref_ids = obs["soma_joinid"].sample(
            min(f["n_reference_cells"], len(obs)),
            random_state=f["seed"],
        )
        reference = cellxgene_census.get_anndata(
            census, "Homo sapiens", obs_coords=ref_ids.tolist()
        )
    finally:
        census.close()

    # Keep only obs columns the pipeline needs; index vars by HGNC symbol
    # (census uses Ensembl feature_id) so config marker panels resolve.
    for a in (organoid, reference):
        if "feature_name" in a.var.columns:
            a.var_names = a.var.pop("feature_name").astype(str)
            a.var_names_make_unique()
            a.var.index.name = "feature_name"
    organoid.obs = organoid.obs[[f["organoid_label_col"]]].copy()
    reference.obs = reference.obs[["cell_type"]].copy()
    return organoid, reference


def main() -> None:
    org_out = sys.argv[1] if len(sys.argv) > 1 else "data/raw/organoid.h5ad"
    ref_out = sys.argv[2] if len(sys.argv) > 2 else "data/raw/reference.h5ad"
    cfg = load_config()
    organoid, reference = fetch(cfg)
    organoid.write_h5ad(org_out)
    reference.write_h5ad(ref_out)
    print(
        f"fetch: organoid {organoid.n_obs} cells -> {org_out}; "
        f"reference {reference.n_obs} cells -> {ref_out}"
    )
    print(
        "set dataset.mode: h5ad and point organoid_h5ad / reference_h5ad "
        "at these files"
    )


if __name__ == "__main__":
    main()
