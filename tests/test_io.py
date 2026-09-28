from pathlib import Path

import pytest

from cytof_qc import io

DATA_DIR = Path(__file__).resolve().parent.parent / "data" / "Levine_13dim"


@pytest.mark.parametrize(
    ("name", "expected"),
    [
        ("Marrow1_Naive CD8+ T_cells.fcs", ("Naive CD8+ T", None)),
        ("Marrow1_NotGated_cells.fcs", ("NotGated", None)),
        (
            "levine32_AML08_normalized_Basophils_H1.fcs",
            ("Basophils", "H1"),
        ),
        (
            "levine32_AML09_normalized_CD4+ T_H2.fcs",
            ("CD4+ T", "H2"),
        ),
        (
            "levine32_AML08_normalized_NotDebrisSinglets_H1.fcs",
            ("NotGated", "H1"),
        ),
        ("unparseable_name.fcs", ("unknown", None)),
    ],
)
def test_parse_filename(name, expected):
    assert io._parse_filename(name) == expected


@pytest.mark.skipif(not DATA_DIR.exists(), reason="run scripts/fetch_data.sh first")
def test_load_fcs_renames_channels_to_markers():
    fcs = sorted(DATA_DIR.glob("*.fcs"))[0]
    df, meta = io.load_fcs(fcs)
    assert "CD45" in df.columns
    assert not df.empty
    assert meta["$CYT"] == "DVSSCIENCES-CYTOF"


@pytest.mark.skipif(not DATA_DIR.exists(), reason="run scripts/fetch_data.sh first")
def test_load_population_dir_labels_from_filename():
    events = io.load_population_dir(DATA_DIR)
    assert "population" in events.columns
    assert "Naive CD8+ T" in set(events["population"])
    assert "NotGated" in set(events["population"])
    assert events["population"].nunique() >= 20
