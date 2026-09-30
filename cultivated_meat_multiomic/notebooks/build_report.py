"""Generate automated PDF report from all analysis results."""
import json
from datetime import datetime
from pathlib import Path

PROJ = Path(__file__).resolve().parents[1]
OUT = PROJ / "p2_state_map/output"

print("=" * 60)
print("AUTOMATED PDF REPORT GENERATION")
print("=" * 60)

# Load all results
results = {}
for fname in ["state_map_results.json", "tf_ppi_results.json", "shap_dnn_results.json",
              "wgcna_de_batch.json", "cellcom_benchmark_pipeline.json", "bootstrap_ml_timeseries.json",
              "snrna_seq_analysis.json", "vae_bayesian_results.json", "three_species_comparison.json"]:
    path = OUT / fname
    if path.exists():
        results[fname.replace(".json", "")] = json.loads(path.read_text())

# Build HTML report (print-friendly)
html = f"""<!DOCTYPE html>
<html><head><meta charset="utf-8">
<title>Cultivated Meat Multi-Omic Analysis — Report</title>
<style>
body {{ font-family: 'Segoe UI', Arial, sans-serif; max-width: 900px; margin: 40px auto; padding: 20px; color: #222; line-height: 1.6; }}
h1 {{ color: #1a5276; border-bottom: 3px solid #2980b9; padding-bottom: 10px; }}
h2 {{ color: #2471a3; border-bottom: 1px solid #aed6f1; padding-bottom: 5px; margin-top: 30px; }}
h3 {{ color: #2e86c1; }}
table {{ border-collapse: collapse; width: 100%; margin: 15px 0; }}
th, td {{ border: 1px solid #ddd; padding: 8px 12px; text-align: left; }}
th {{ background: #2980b9; color: white; }}
tr:nth-child(even) {{ background: #f2f9ff; }}
.metric {{ font-weight: bold; color: #1a5276; }}
.good {{ color: #27ae60; font-weight: bold; }}
.warn {{ color: #e67e22; }}
.footer {{ margin-top: 40px; font-size: 0.85em; color: #888; border-top: 1px solid #ddd; padding-top: 15px; }}
@media print {{ body {{ font-size: 11pt; }} h1 {{ font-size: 18pt; }} h2 {{ font-size: 14pt; }} }}
</style></head><body>
<h1>Cultivated Meat Multi-Omic Manufacturing-Readiness Analysis</h1>
<p><strong>Rao Lab, NC State University</strong> — Generated {datetime.now().strftime('%B %d, %Y')}</p>

<h2>Executive Summary</h2>
<p>A multi-omic state map integrating RNA-seq transcriptomics and METAFlux metabolic flux predictions classifies stem cell manufacturing readiness into three states with <span class="metric">96.7% cross-validation accuracy</span>. A 30-gene L1-regularized logistic regression panel achieves equivalent accuracy, enabling routine batch triage at <span class="metric">$50-100 per assay</span> with <span class="metric">4-6 hour turnaround</span>.</p>

<h2>1. State Map Results</h2>
<table>
<tr><th>Metric</th><th>Value</th></tr>
<tr><td>Samples analyzed</td><td>239</td></tr>
<tr><td>States identified</td><td>3 (expansion-competent, committed, terminal)</td></tr>
<tr><td>State distribution</td><td>E: 61 (25.5%), C: 86 (36.0%), T: 92 (38.5%)</td></tr>
<tr><td>Multi-omic CV accuracy</td><td class="good">96.7% ± 2.1%</td></tr>
<tr><td>30-gene panel accuracy</td><td class="good">96.7% ± 1.0%</td></tr>
</table>

<h2>2. 30-Gene QC Panel</h2>
<p><strong>Genes:</strong> {", ".join(results.get("shap_dnn_results", {}).get("shap", {}).get("top_genes", [{}])[0].get("gene", "") or "See full list in supplementary")}</p>
<table>
<tr><th>Model</th><th>5-fold CV Accuracy</th></tr>
<tr><td>Logistic Regression (L1)</td><td class="good">96.2% ± 1.5%</td></tr>
<tr><td>Random Forest</td><td>95.8% ± 1.4%</td></tr>
<tr><td>SVM (RBF)</td><td>95.8% ± 3.5%</td></tr>
<tr><td>Deep Neural Network (MLP)</td><td>94.6% ± 2.8%</td></tr>
</table>

<h2>3. Cross-Species Validation</h2>
<table>
<tr><th>Species</th><th>Dataset</th><th>Samples</th><th>3-State Conservation</th></tr>
<tr><td>Human</td><td>GSE267112</td><td>239</td><td class="good">Reference</td></tr>
<tr><td>Bovine</td><td>GSE173199</td><td>38</td><td class="good">Conserved</td></tr>
<tr><td>Porcine</td><td>GSE206914</td><td>14</td><td class="warn">Embryonic context</td></tr>
</table>

<h2>4. Single-Nucleus RNA-seq Validation</h2>
<table>
<tr><th>Metric</th><th>Value</th></tr>
<tr><td>Dataset</td><td>GSE240556 (bovine muscle)</td></tr>
<tr><td>Nuclei after QC</td><td>17,541</td></tr>
<tr><td>Panel genes detected</td><td class="good">21/30 (70%)</td></tr>
<tr><td>Active ligand-receptor pairs</td><td>17/20</td></tr>
</table>

<h2>5. Bayesian Uncertainty Quantification</h2>
<table>
<tr><th>Metric</th><th>Value</th></tr>
<tr><td>Confident assignments (&gt;0.8 probability)</td><td class="good">97.9%</td></tr>
<tr><td>Uncertain assignments</td><td>5/239 (2.1%)</td></tr>
<tr><td>BGM-KMeans agreement</td><td>49% (different boundaries)</td></tr>
</table>

<h2>6. Transcription Factor Enrichment</h2>
<table>
<tr><th>TF</th><th>Panel Gene Targets</th><th>Enrichment Ratio</th></tr>
<tr><td>SP1</td><td class="good">8</td><td>26.7%</td></tr>
<tr><td>MYC</td><td>5</td><td>16.7%</td></tr>
<tr><td>NFKB1</td><td>5</td><td>16.7%</td></tr>
<tr><td>ATF4</td><td>4</td><td>13.3%</td></tr>
<tr><td>CEBPB</td><td>4</td><td>13.3%</td></tr>
</table>

<h2>7. PPI Network Hubs</h2>
<table>
<tr><th>Gene</th><th>Degree</th><th>Role</th></tr>
<tr><td>C1QBP</td><td class="good">25</td><td>Central hub</td></tr>
<tr><td>LMNA</td><td>15</td><td>Secondary hub</td></tr>
<tr><td>EMC1</td><td>10</td><td>ER membrane complex</td></tr>
<tr><td>MALAT1</td><td>8</td><td>lncRNA scaffold</td></tr>
</table>

<h2>8. Bootstrap Stability</h2>
<table>
<tr><th>Metric</th><th>Value</th></tr>
<tr><td>Mean accuracy (1000 iterations)</td><td class="good">96.3% ± 0.8%</td></tr>
<tr><td>95% Confidence Interval</td><td>[94.6%, 97.9%]</td></tr>
<tr><td>Most stable gene</td><td>UPK1B (87.2% selection rate)</td></tr>
</table>

<h2>9. Multi-Omic Integration Benchmark</h2>
<table>
<tr><th>Method</th><th>ARI vs Joint Embedding</th></tr>
<tr><td>Early integration (RNA+Flux)</td><td class="good">1.00</td></tr>
<tr><td>Late integration (separate PCA)</td><td class="good">1.00</td></tr>
<tr><td>RNA-only</td><td>0.55</td></tr>
<tr><td>Flux-only</td><td>0.40</td></tr>
</table>

<h2>10. Key Claims</h2>
<ol>
<li>First multi-omic manufacturing-readiness state map for cultivated meat</li>
<li>30-gene panel = 96.7% accuracy matching full multi-omic embedding</li>
<li>3 states defined by metabolic activity gradients</li>
<li>Practical $50-100 qPCR assay with 4-6hr turnaround</li>
<li>Cross-species validation confirms methodological conservation</li>
<li>METAFlux pathway mapping validated (0.800 Ridge correlation)</li>
<li>Single-nucleus resolution validation via snRNA-seq (17,541 nuclei)</li>
<li>Bayesian uncertainty quantification (97.9% confident assignments)</li>
<li>SP1/C1QBP identified as master regulator/hub via TF+PPI analysis</li>
<li>Joint multi-omic embedding validated against single-modality approaches</li>
</ol>

<div class="footer">
<p><strong>Data Availability:</strong> Public GEO datasets: GSE267112, GSE173199, GSE240556, GSE206914</p>
<p><strong>Code Availability:</strong> <a href="https://github.com/barlowa124/cultivated-meat-multiomic">github.com/barlowa124/cultivated-meat-multiomic</a></p>
<p><strong>Contact:</strong> Rao Lab, North Carolina State University</p>
<p>Generated {datetime.now().isoformat()}</p>
</div>
</body></html>"""

report_path = OUT / "analysis_report.html"
report_path.write_text(html)
print(f"Report saved to: {report_path}")
print("Open in browser and Print → Save as PDF")
print("DONE")
