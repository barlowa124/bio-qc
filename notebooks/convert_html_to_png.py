"""Convert interactive Plotly HTML figures to static PNG using Playwright.

More reliable than kaleido on Windows. Opens each HTML file in headless
Chromium and takes a full-page screenshot at high resolution.
"""
import warnings
from pathlib import Path

from playwright.sync_api import sync_playwright

warnings.filterwarnings("ignore")

PROJ = Path(__file__).resolve().parents[1]
FIGS = PROJ / "docs/figures"

HTML_FIGS = [
    "fig1_state_map.html",
    "fig2_qc_performance.html",
    "fig3_shap_importance.html",
    "fig4_cross_species.html",
    "fig5_noise_robustness.html",
    "fig6_drug_predictions.html",
    "fig7_tea_sensitivity.html",
    "fig9_batch_correction.html",
    "fig10_pathway_enrichment.html",
    "fig11_cross_platform.html",
    "fig11b_platform_heatmap.html",
]

print("=" * 60)
print("STATIC FIGURES via Playwright")
print("=" * 60)

with sync_playwright() as p:
    browser = p.chromium.launch(headless=True)
    context = browser.new_context(
        viewport={"width": 1400, "height": 1000},
        device_scale_factor=2,
    )

    for html_file in HTML_FIGS:
        html_path = FIGS / html_file
        png_name = html_file.replace(".html", ".png")

        if not html_path.exists():
            print(f"  SKIP {html_file} (not found)")
            continue

        print(f"  Rendering {html_file} -> {png_name}")

        page = context.new_page()
        page.goto(f"file:///{html_path.as_posix()}", wait_until="networkidle")
        page.wait_for_timeout(1000)  # let Plotly finish rendering

        page.screenshot(
            path=str(FIGS / png_name),
            full_page=True,
        )
        page.close()

    browser.close()

print("\n--- Done ---")
for png in sorted(FIGS.glob("fig*.png")):
    print(f"  {png.name}")
print(f"Output: {FIGS}")
print("DONE")
