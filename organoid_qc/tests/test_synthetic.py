import anndata as ad
import yaml

from organoid_qc.data.synthetic import make_demo


def _cfg():
    with open("config/config.yaml") as f:
        return yaml.safe_load(f)


def test_demo_shapes_and_types():
    organoid, reference = make_demo(_cfg())
    assert isinstance(organoid, ad.AnnData)
    assert reference.n_obs == 900
    assert organoid.n_obs == 600
    assert set(reference.obs["cell_type"]) == {"enterocyte", "goblet", "stem"}
    # Ground truth: one type missing, one aberrant population present.
    assert organoid.uns["demo_missing_type"] == "stem"
    assert "stem" not in set(organoid.obs["true_population"])
    assert "unannotated" in set(organoid.obs["true_population"])


def test_demo_is_deterministic():
    o1, r1 = make_demo(_cfg())
    o2, r2 = make_demo(_cfg())
    assert (o1.X == o2.X).all()
    assert (r1.X == r2.X).all()
