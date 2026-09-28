# Agent notes

## Testing

- Primary suite: `.venv/bin/python -m pytest -q` (Python 3.9 with pytest, flask, sklearn. No plotly or streamlit, so figure and dashboard tests skip there.)
- Plotly/Streamlit coverage: `.venv313` (Python 3.13) has plotly, pandas, streamlit AppTest, kaleido, but no pytest or flask. Symlink the pure-Python pytest packages (pytest, _pytest, pluggy, iniconfig, packaging, pygments, py.py) from `.venv/lib/python3.9/site-packages` into a staging dir such as `/tmp/pytest313`, then run `PYTHONPATH=/tmp/pytest313 .venv313/bin/python -m pytest tests/test_provenance.py -q`. Do not point PYTHONPATH at the whole 3.9 site-packages: its numpy C-extensions fail to import under 3.13.

## Warning: fixture and generator scripts overwrite committed artifacts

- `ci_fixtures.py` writes `api/`, `p2_state_map/output/`, and `p3_qc_panel/output/` relative to its own directory. Never run it at the repo root. Copy it into a temp directory and run it there.
- `notebooks/generate_interactive_figures.py`, `generate_static_figures.py`, and `generate_svg_figures.py` write into `docs/figures/` relative to the repo root. Run copies staged in a temp directory, or accept that committed figures will be regenerated.
