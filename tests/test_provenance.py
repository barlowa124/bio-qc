import ast
import importlib.util
import json
import runpy
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

PROJ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJ / "notebooks"))

import artifact_provenance as ap  # noqa: E402

GEN_SCRIPTS = [
    "generate_interactive_figures.py",
    "generate_static_figures.py",
    "generate_svg_figures.py",
]


@pytest.mark.parametrize("source", ["ci_fixture_synthetic", "synthetic_formula", "simulated_subset"])
def test_synthetic_sources_flagged(source):
    artifact = {"data_source": source}
    assert ap.is_synthetic(artifact) is True
    label = ap.provenance_label(artifact)
    assert "Synthetic data" in label
    assert source in label


@pytest.mark.parametrize(
    "artifact",
    [{}, {"data_source": ""}, {"data_source": "   "}, {"data_source": 7},
     ["not", "a", "dict"], None],
)
def test_missing_source_is_unverified(artifact):
    assert ap.data_source(artifact) == "unverified"
    assert ap.is_synthetic(artifact) is False
    label = ap.provenance_label(artifact)
    assert "Unverified" in label
    for word in ("measured", "real", "Recorded"):
        assert word not in label


def test_recorded_table_source():
    artifact = {"data_source": "Table_S2_ml_benchmark.csv"}
    assert ap.is_synthetic(artifact) is False
    assert ap.provenance_label(artifact) == "Recorded source: Table_S2_ml_benchmark.csv"


def test_annotate_figure_marks_title_and_meta():
    go = pytest.importorskip("plotly.graph_objects")
    fig = go.Figure()
    fig.update_layout(title_text="Base Title", meta={"custom_key": "keep"})
    out = ap.annotate_figure(fig, {"data_source": "ci_fixture_synthetic"})
    assert out is fig
    ap.annotate_figure(fig, {"data_source": "ci_fixture_synthetic"})
    title = fig.layout.title.text
    assert title.count("Synthetic data (ci_fixture_synthetic)") == 1
    assert fig.layout.meta["custom_key"] == "keep"
    assert fig.layout.meta["data_source"] == "ci_fixture_synthetic"
    assert fig.layout.meta["synthetic"] is True


def test_annotate_figure_unverified_label():
    go = pytest.importorskip("plotly.graph_objects")
    fig = go.Figure()
    fig.update_layout(title_text="T")
    ap.annotate_figure(fig, {})
    assert "Unverified provenance" in fig.layout.title.text
    assert fig.layout.meta["data_source"] == "unverified"
    assert fig.layout.meta["synthetic"] is False


def _generate_ci_fixture_repo(root):
    shutil.copy(PROJ / "ci_fixtures.py", root / "ci_fixtures.py")
    api_dir = root / "api"
    api_dir.mkdir(parents=True, exist_ok=True)
    if not (api_dir / "app.py").exists():
        shutil.copy(PROJ / "api" / "app.py", api_dir / "app.py")
    proc = subprocess.run([sys.executable, str(root / "ci_fixtures.py")],
                          cwd=root, capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr


def test_ci_fixtures_mark_every_artifact(tmp_path):
    _generate_ci_fixture_repo(tmp_path)
    jsons = sorted((tmp_path / "p2_state_map" / "output").glob("*.json"))
    jsons += sorted((tmp_path / "p3_qc_panel" / "output").glob("*.json"))
    jsons.append(tmp_path / "api" / "model_metadata.json")
    assert len(jsons) >= 10
    for path in jsons:
        artifact = json.loads(path.read_text())
        assert artifact.get("data_source") == "ci_fixture_synthetic", path.name
    grant = (tmp_path / "p2_state_map" / "output" / "grant_proposal_draft.md").read_text()
    assert grant.startswith("Synthetic CI fixture.")


def test_production_api_with_fixture_data_reports_provenance(tmp_path):
    pytest.importorskip("flask")
    _generate_ci_fixture_repo(tmp_path)
    spec = importlib.util.spec_from_file_location(
        "fixture_api_app", tmp_path / "api" / "app.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    client = module.app.test_client()
    index = client.get("/").get_json()
    assert index["data_source"] == "ci_fixture_synthetic"
    health = client.get("/health").get_json()
    assert health["model_accuracy"] is None
    assert health["data_source"] == "ci_fixture_synthetic"
    response = client.post("/predict", json={"expression": [1.0] * 30})
    assert response.status_code == 200
    assert response.get_json()["data_source"] == "ci_fixture_synthetic"


def test_fixture_generator_preserves_api_source(tmp_path):
    api_dir = tmp_path / "api"
    api_dir.mkdir()
    sentinel = api_dir / "app.py"
    sentinel.write_bytes(b"# unique sentinel 7f3a9c: fixture must not overwrite\n")
    _generate_ci_fixture_repo(tmp_path)
    assert sentinel.read_bytes() == b"# unique sentinel 7f3a9c: fixture must not overwrite\n"


def test_fixture_generator_does_not_create_api_source(tmp_path):
    shutil.copy(PROJ / "ci_fixtures.py", tmp_path / "ci_fixtures.py")
    proc = subprocess.run([sys.executable, str(tmp_path / "ci_fixtures.py")],
                          cwd=tmp_path, capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    assert not (tmp_path / "api" / "app.py").exists()


def _stage_notebooks(root):
    nb = root / "notebooks"
    nb.mkdir(parents=True, exist_ok=True)
    for name in GEN_SCRIPTS + ["artifact_provenance.py"]:
        shutil.copy(PROJ / "notebooks" / name, nb / name)
    (root / "p2_state_map" / "output").mkdir(parents=True, exist_ok=True)
    (root / "docs" / "figures").mkdir(parents=True, exist_ok=True)
    (root / "docs" / "supplementary").mkdir(parents=True, exist_ok=True)


@pytest.fixture(scope="module")
def fixture_repo(tmp_path_factory):
    root = tmp_path_factory.mktemp("fixture_repo")
    _generate_ci_fixture_repo(root)
    _stage_notebooks(root)
    return root


@pytest.fixture
def captured_figs(monkeypatch):
    pytest.importorskip("plotly.graph_objects")
    from plotly.basedatatypes import BaseFigure

    captured = {}

    def capture(self, file, *args, **kwargs):
        captured[Path(file).name] = self.to_dict()

    monkeypatch.setattr(BaseFigure, "write_html", capture)
    monkeypatch.setattr(BaseFigure, "write_image", capture)
    return captured


def _fig(captured, stem):
    for name, fig_dict in captured.items():
        if Path(name).stem == stem:
            return fig_dict
    return None


def _run_generator(root, script, monkeypatch):
    monkeypatch.syspath_prepend(str(root / "notebooks"))
    runpy.run_path(str(root / "notebooks" / script), run_name="__main__")


FIXTURE_EXPECTED = {
    "generate_interactive_figures.py": {
        "fig1_state_map": ("ci_fixture_synthetic", True),
        "fig2_qc_performance": ("unverified", False),
        "fig3_shap_importance": ("ci_fixture_synthetic", True),
        "fig4_cross_species": ("ci_fixture_synthetic", True),
        "fig5_noise_robustness": ("ci_fixture_synthetic", True),
        "fig6_drug_predictions": ("ci_fixture_synthetic", True),
        "fig7_tea_sensitivity": ("synthetic_formula", True),
        "fig8_ppi_network": ("ci_fixture_synthetic", True),
    },
    "generate_static_figures.py": {
        "fig1_state_map": ("ci_fixture_synthetic", True),
        "fig3_shap_importance": ("ci_fixture_synthetic", True),
        "fig4_cross_species": ("ci_fixture_synthetic", True),
        "fig5_noise_robustness": ("ci_fixture_synthetic", True),
        "fig6_drug_predictions": ("ci_fixture_synthetic", True),
        "fig7_tea_sensitivity": ("synthetic_formula", True),
        "fig8_ppi_network": ("ci_fixture_synthetic", True),
    },
    "generate_svg_figures.py": {
        "fig1_state_map": ("ci_fixture_synthetic", True),
        "fig3_shap_importance": ("ci_fixture_synthetic", True),
        "fig4_cross_species": ("ci_fixture_synthetic", True),
        "fig5_noise_robustness": ("ci_fixture_synthetic", True),
        "fig6_drug_predictions": ("ci_fixture_synthetic", True),
        "fig7_tea_sensitivity": ("synthetic_formula", True),
    },
}


@pytest.mark.parametrize("script", sorted(FIXTURE_EXPECTED))
def test_generators_label_fixture_sources(fixture_repo, captured_figs, monkeypatch, script):
    _run_generator(fixture_repo, script, monkeypatch)
    for stem, (source, synthetic) in FIXTURE_EXPECTED[script].items():
        fig_dict = _fig(captured_figs, stem)
        assert fig_dict is not None, f"{script} did not emit {stem}"
        meta = fig_dict["layout"]["meta"]
        assert meta["data_source"] == source, stem
        assert meta["synthetic"] is synthetic, stem
        title = fig_dict["layout"]["title"]["text"]
        if synthetic:
            assert "Synthetic data" in title, stem
        else:
            assert "Unverified provenance" in title, stem
    go = pytest.importorskip("plotly.graph_objects")
    html = go.Figure(_fig(captured_figs, "fig1_state_map")).to_html(include_plotlyjs="cdn")
    assert "Synthetic data" in html
    assert "ci_fixture_synthetic" in html


EMPTY_EXPECTED = {
    "generate_interactive_figures.py": {
        "emitted": {
            "fig2_qc_performance": "unverified",
            "fig4_cross_species": "synthetic_fallback",
            "fig7_tea_sensitivity": "synthetic_formula",
        },
        "absent": ["fig1_state_map", "fig3_shap_importance", "fig5_noise_robustness",
                   "fig6_drug_predictions", "fig8_ppi_network"],
    },
    "generate_static_figures.py": {
        "emitted": {
            "fig4_cross_species": "synthetic_fallback",
            "fig5_noise_robustness": "synthetic_formula",
            "fig7_tea_sensitivity": "synthetic_formula",
        },
        "absent": ["fig1_state_map", "fig2_qc_performance", "fig3_shap_importance",
                   "fig6_drug_predictions", "fig8_ppi_network"],
    },
    "generate_svg_figures.py": {
        "emitted": {
            "fig4_cross_species": "synthetic_fallback",
            "fig5_noise_robustness": "synthetic_formula",
            "fig7_tea_sensitivity": "synthetic_formula",
        },
        "absent": ["fig1_state_map", "fig2_qc_performance", "fig3_shap_importance",
                   "fig6_drug_predictions"],
    },
}


@pytest.mark.parametrize("script", sorted(EMPTY_EXPECTED))
def test_generators_label_fallbacks(tmp_path, captured_figs, monkeypatch, script):
    _stage_notebooks(tmp_path)
    _run_generator(tmp_path, script, monkeypatch)
    spec = EMPTY_EXPECTED[script]
    for stem, source in spec["emitted"].items():
        fig_dict = _fig(captured_figs, stem)
        assert fig_dict is not None, f"{script} did not emit {stem}"
        meta = fig_dict["layout"]["meta"]
        assert meta["data_source"] == source, stem
        title = fig_dict["layout"]["title"]["text"]
        if source == "unverified":
            assert meta["synthetic"] is False
            assert "Unverified provenance" in title
        else:
            assert meta["synthetic"] is True
            assert "Synthetic data" in title
    for stem in spec["absent"]:
        assert _fig(captured_figs, stem) is None, f"{script} should not emit {stem}"


@pytest.mark.parametrize("script", GEN_SCRIPTS)
def test_unmarked_artifacts_label_unverified(tmp_path, captured_figs, monkeypatch, script):
    _stage_notebooks(tmp_path)
    (tmp_path / "p2_state_map" / "output" / "state_map_results.json").write_text(json.dumps({
        "n_samples": 6,
        "states": {"expansion_competent": 2, "committed": 2, "terminal": 2},
        "accuracy": 0.9,
    }))
    _run_generator(tmp_path, script, monkeypatch)
    fig_dict = _fig(captured_figs, "fig1_state_map")
    assert fig_dict is not None, f"{script} did not emit fig1"
    meta = fig_dict["layout"]["meta"]
    assert meta["data_source"] == "unverified"
    assert meta["synthetic"] is False
    title = fig_dict["layout"]["title"]["text"]
    assert "Unverified provenance" in title
    assert "Recorded source" not in title
    assert "Synthetic" not in title


TABLES = [
    "Table_S2_ml_benchmark.csv",
    "Table_S4_cross_platform.csv",
    "Table_S5_batch_correction.csv",
    "Table_S6_pathway_enrichment.csv",
]


def _write_synthetic_domain_jsons(out, concordance):
    (out / "batch_correction_benchmark.json").write_text(json.dumps({
        "data_source": "ci_fixture_synthetic",
        "methods": {"raw": {"accuracy": 0.9, "batch_mixing": 0.4}},
    }))
    (out / "pathway_enrichment.json").write_text(json.dumps({
        "data_source": "ci_fixture_synthetic",
        "go_bp": {"top_hits": [{"term": "lipid_metabolism", "p_value": 0.01,
                                "p_value_adj": 0.02, "overlap_count": 2, "term_size": 10}]},
    }))
    (out / "cross_platform_validation.json").write_text(json.dumps({
        "data_source": "ci_fixture_synthetic",
        "expression_concordance": concordance,
        "gene_bias_summary": {"LMNA": {"rna_seq_mean": 1.0, "qPCR_mean": 1.0,
                                       "nanostring_mean": 1.0}},
    }))


@pytest.fixture
def table_repo(tmp_path):
    _stage_notebooks(tmp_path)
    supp = tmp_path / "docs" / "supplementary"
    for name in TABLES:
        shutil.copy(PROJ / "docs" / "supplementary" / name, supp / name)
    out = tmp_path / "p2_state_map" / "output"
    _write_synthetic_domain_jsons(out, concordance=[
        {"platform_a": "rna_seq", "platform_b": "qPCR", "mean_gene_correlation": 0.99},
    ])
    return tmp_path


TABLE_EXPECTED = {
    "generate_interactive_figures.py": {
        "fig9_batch_correction": "Table_S5_batch_correction.csv",
        "fig10_pathway_enrichment": "Table_S6_pathway_enrichment.csv",
        "fig11_cross_platform": "Table_S4_cross_platform.csv",
    },
    "generate_static_figures.py": {
        "fig2_qc_performance": "Table_S2_ml_benchmark.csv",
        "fig9_batch_correction": "Table_S5_batch_correction.csv",
        "fig10_pathway_enrichment": "Table_S6_pathway_enrichment.csv",
        "fig11_cross_platform": "Table_S4_cross_platform.csv",
        "fig11b_platform_heatmap": "Table_S4_cross_platform.csv",
    },
    "generate_svg_figures.py": {
        "fig2_qc_performance": "Table_S2_ml_benchmark.csv",
        "fig9_batch_correction": "Table_S5_batch_correction.csv",
        "fig10_pathway_enrichment": "Table_S6_pathway_enrichment.csv",
        "fig11_cross_platform": "Table_S4_cross_platform.csv",
        "fig11b_platform_heatmap": "Table_S4_cross_platform.csv",
    },
}


@pytest.mark.parametrize("script", sorted(TABLE_EXPECTED))
def test_table_preferred_over_synthetic_json(table_repo, captured_figs, monkeypatch, script):
    _run_generator(table_repo, script, monkeypatch)
    for stem, source in TABLE_EXPECTED[script].items():
        fig_dict = _fig(captured_figs, stem)
        assert fig_dict is not None, f"{script} did not emit {stem}"
        meta = fig_dict["layout"]["meta"]
        assert meta["data_source"] == source, stem
        assert meta["synthetic"] is False, stem
        title = fig_dict["layout"]["title"]["text"]
        assert f"Recorded source: {source}" in title
        assert "ci_fixture_synthetic" not in title
    if script == "generate_interactive_figures.py":
        fig2 = _fig(captured_figs, "fig2_qc_performance")
        assert fig2["layout"]["meta"]["data_source"] == "unverified"
        assert _fig(captured_figs, "fig11b_platform_heatmap") is None


TEA_FIXTURE = {
    "data_source": "ci_fixture_synthetic",
    "sensitivity": [{"panel_cost_per_batch": 50, "cost_per_sample": 1},
                    {"panel_cost_per_batch": 100, "cost_per_sample": 2}],
}

FULL_CONCORDANCE = [
    {"platform_a": "rna_seq", "platform_b": "qPCR", "mean_gene_correlation": 0.99},
    {"platform_a": "rna_seq", "platform_b": "nanostring", "mean_gene_correlation": 0.98},
    {"platform_a": "qPCR", "platform_b": "nanostring", "mean_gene_correlation": 0.97},
]


@pytest.mark.parametrize("script", ["generate_interactive_figures.py", "generate_static_figures.py"])
def test_populated_synthetic_json_labeled(tmp_path, captured_figs, monkeypatch, script):
    _stage_notebooks(tmp_path)
    out = tmp_path / "p2_state_map" / "output"
    _write_synthetic_domain_jsons(out, concordance=FULL_CONCORDANCE)
    (out / "techno_economic_analysis.json").write_text(json.dumps(TEA_FIXTURE))
    _run_generator(tmp_path, script, monkeypatch)
    for stem in ("fig7_tea_sensitivity", "fig9_batch_correction", "fig10_pathway_enrichment",
                 "fig11_cross_platform", "fig11b_platform_heatmap"):
        fig_dict = _fig(captured_figs, stem)
        assert fig_dict is not None, f"{script} did not emit {stem}"
        meta = fig_dict["layout"]["meta"]
        assert meta["data_source"] == "ci_fixture_synthetic", stem
        assert meta["synthetic"] is True, stem
        title = fig_dict["layout"]["title"]["text"]
        assert "Synthetic data" in title, stem
        assert "Table S" not in title, stem
    titles = {stem: _fig(captured_figs, stem)["layout"]["title"]["text"]
              for stem in ("fig9_batch_correction", "fig10_pathway_enrichment",
                           "fig11_cross_platform", "fig11b_platform_heatmap")}
    assert "pathway_enrichment.json" in titles["fig10_pathway_enrichment"]
    if script == "generate_static_figures.py":
        assert "batch_correction_benchmark.json" in titles["fig9_batch_correction"]
        assert "cross_platform_validation.json" in titles["fig11_cross_platform"]
        assert "cross_platform_validation.json" in titles["fig11b_platform_heatmap"]


def test_svg_populated_tea_labeled(tmp_path, captured_figs, monkeypatch):
    _stage_notebooks(tmp_path)
    (tmp_path / "p2_state_map" / "output" / "techno_economic_analysis.json").write_text(
        json.dumps(TEA_FIXTURE))
    _run_generator(tmp_path, "generate_svg_figures.py", monkeypatch)
    fig_dict = _fig(captured_figs, "fig7_tea_sensitivity")
    assert fig_dict is not None
    meta = fig_dict["layout"]["meta"]
    assert meta["data_source"] == "ci_fixture_synthetic"
    assert meta["synthetic"] is True
    assert "Synthetic data" in fig_dict["layout"]["title"]["text"]


def _stage_dashboard_repo(root):
    _generate_ci_fixture_repo(root)
    nb = root / "notebooks"
    nb.mkdir(parents=True, exist_ok=True)
    shutil.copy(PROJ / "notebooks" / "artifact_provenance.py", nb / "artifact_provenance.py")
    dash = root / "dashboard"
    dash.mkdir(parents=True, exist_ok=True)
    shutil.copy(PROJ / "dashboard" / "app.py", dash / "app.py")
    return dash / "app.py"


def test_dashboard_warns_and_annotates(tmp_path):
    pytest.importorskip("streamlit")
    from streamlit.testing.v1 import AppTest

    app_path = _stage_dashboard_repo(tmp_path)
    at = AppTest.from_file(str(app_path), default_timeout=120)
    at.run()
    assert not at.exception
    warnings = [w.value for w in at.warning]
    assert any("model_metadata.json" in w for w in warnings)
    assert any("shap_dnn_results.json" in w for w in warnings)
    assert any("Synthetic data" in w for w in warnings)

    at.button[0].click().run()
    assert not at.exception
    charts = at.get("plotly_chart")
    assert charts
    spec = str(charts[0].proto)
    assert "ci_fixture_synthetic" in spec
    assert "Synthetic data" in spec

    at.sidebar.radio[0].set_value("Performance").run()
    assert not at.exception
    perf_warnings = [w.value for w in at.warning]
    assert any("bootstrap_ml_timeseries.json" in w for w in perf_warnings)
    assert all("Synthetic data" in w for w in perf_warnings)
    charts = at.get("plotly_chart")
    assert charts
    assert "ci_fixture_synthetic" in str(charts[0].proto)

    at.sidebar.radio[0].set_value("Explainability").run()
    assert not at.exception
    metrics = [m.value for m in at.metric]
    assert "LMNA" in metrics

    at.sidebar.radio[0].set_value("Literature & Drugs").run()
    assert not at.exception
    lit_warnings = [w.value for w in at.warning]
    assert any("literature_drug_panel_noise.json" in w for w in lit_warnings)
    assert all("Synthetic data" in w for w in lit_warnings)


def test_dashboard_unverified_model_metadata(tmp_path):
    pytest.importorskip("streamlit")
    import streamlit
    from streamlit.testing.v1 import AppTest

    app_path = _stage_dashboard_repo(tmp_path)
    meta_path = tmp_path / "api" / "model_metadata.json"
    meta = json.loads(meta_path.read_text())
    meta.pop("data_source", None)
    meta_path.write_text(json.dumps(meta))
    streamlit.cache_resource.clear()
    streamlit.cache_data.clear()
    at = AppTest.from_file(str(app_path), default_timeout=120)
    at.run()
    assert not at.exception
    warnings = [w.value for w in at.warning]
    assert any("model_metadata.json: Unverified provenance" in w for w in warnings)
    assert not any("model_metadata.json" in w and "Recorded source" in w for w in warnings)


def test_dashboard_batch_export_marks_source():
    tree = ast.parse((PROJ / "dashboard" / "app.py").read_text())
    hits = [
        node for node in ast.walk(tree)
        if isinstance(node, ast.Assign)
        and any(
            isinstance(t, ast.Subscript)
            and isinstance(t.slice, ast.Constant)
            and t.slice.value == "model_data_source"
            for t in node.targets
        )
        and isinstance(node.value, ast.Call)
        and getattr(node.value.func, "id", "") == "data_source"
    ]
    assert hits
