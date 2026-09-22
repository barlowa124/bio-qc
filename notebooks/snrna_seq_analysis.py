"""GSE240556 bovine snRNA-seq analysis with Scanpy."""
import gzip
import io
import json
import tarfile
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
PROJ = Path(__file__).resolve().parents[1]
CROSS = PROJ / "cross_species_validation"
OUT = PROJ / "p2_state_map/output"
OUT.mkdir(exist_ok=True)

print("=" * 60)
print("GSE240556 BOVINE snRNA-seq ANALYSIS")
print("=" * 60)

gse_dir = CROSS / "GSE240556"
suppl_files = list(gse_dir.glob("suppl_*.txt"))

print(f"\nFound {len(suppl_files)} supplementary files")

# These are tar archives containing 10x Genomics output:
# barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz
# Let's try to extract them

extracted = {}
for sf in suppl_files:
    print(f"\nProcessing {sf.name}...")
    try:
        data = sf.read_bytes()
        extracted_any = False
        
        # Try as raw tar (no gzip)
        try:
            with tarfile.open(fileobj=io.BytesIO(data)) as tar:
                members = tar.getmembers()
                print(f"  Raw tar with {len(members)} members:")
                for m in members:
                    print(f"    {m.name} ({m.size} bytes)")
                    f = tar.extractfile(m)
                    if f:
                        content = f.read()
                        extracted[m.name] = content
                        extracted_any = True
                        print(f"      -> extracted {len(content)} bytes")
        except (tarfile.TarError, OSError, EOFError):
            pass
        
        # Try as gzipped tar
        if not extracted_any:
            try:
                with gzip.GzipFile(fileobj=io.BytesIO(data)) as gz:
                    tar_bytes = gz.read()
                with tarfile.open(fileobj=io.BytesIO(tar_bytes)) as tar:
                    members = tar.getmembers()
                    print(f"  Gzipped tar with {len(members)} members:")
                    for m in members:
                        print(f"    {m.name} ({m.size} bytes)")
                        f = tar.extractfile(m)
                        if f:
                            content = f.read()
                            extracted[m.name] = content
                            extracted_any = True
                            print(f"      -> extracted {len(content)} bytes")
            except (tarfile.TarError, gzip.BadGzipFile, OSError, EOFError):
                pass
        
        # Try as direct gzip
        if not extracted_any:
            try:
                content = gzip.decompress(data)
                print(f"  Direct gzip: {len(content)} bytes")
                text = content.decode("utf-8", errors="ignore")[:500]
                print(f"  Preview: {text[:200]}")
                extracted_any = True
            except (gzip.BadGzipFile, OSError, EOFError):
                pass
        
        if not extracted_any:
            print("  Not an archive - likely FTP listing")
            text = data.decode("utf-8", errors="ignore")[:300]
            print(f"  Preview: {text[:200]}")
    except Exception as e:
        print(f"  Error: {e}")

# ── If we got 10x data, process with anndata ──
if extracted:
    print("\n─── Processing extracted 10x data ───")
    try:
        import anndata
        import scipy.io
        import scipy.sparse

        # Find matrix, barcodes, features
        mtx_key = [k for k in extracted if 'mtx' in k.lower()]
        barcode_key = [k for k in extracted if 'barcode' in k.lower()]
        feature_key = [k for k in extracted if 'feature' in k.lower() or 'gene' in k.lower()]

        if mtx_key and barcode_key and feature_key:
            # Decompress and parse
            mtx_data = gzip.decompress(extracted[mtx_key[0]])
            mtx = scipy.io.mmread(io.BytesIO(mtx_data)).tocsr()

            bc_data = gzip.decompress(extracted[barcode_key[0]])
            barcodes = [l.decode().strip() for l in bc_data.split(b'\n') if l.strip()]

            ft_data = gzip.decompress(extracted[feature_key[0]])
            features = [l.decode().strip().split('\t') for l in ft_data.split(b'\n') if l.strip()]

            print(f"  Matrix: {mtx.shape}")
            print(f"  Barcodes: {len(barcodes)}")
            print(f"  Features: {len(features)}")

            # Create AnnData
            adata = anndata.AnnData(
                X=mtx.T.tocsr(),
                obs=pd.DataFrame(index=barcodes),
                var=pd.DataFrame(features, columns=['gene_id','gene_name','feature_type']).set_index('gene_id')
            )
            adata.var_names = [f[1] for f in features]

            # Basic QC
            adata.obs['n_genes'] = np.array((adata.X > 0).sum(axis=1)).flatten()
            adata.obs['n_counts'] = np.array(adata.X.sum(axis=1)).flatten()
            adata.var['n_cells'] = np.array((adata.X > 0).sum(axis=0)).flatten()

            print("\n  QC stats:")
            print(f"    Median genes/cell: {np.median(adata.obs['n_genes']):.0f}")
            print(f"    Median counts/cell: {np.median(adata.obs['n_counts']):.0f}")
            print(f"    Total cells: {adata.n_obs}")

            # Filter
            adata = adata[adata.obs['n_genes'] > 200, :]
            adata = adata[adata.obs['n_genes'] < 5000, :]
            print(f"    After filtering: {adata.n_obs} cells")

            # Check 30-gene panel overlap
            qc_panel_path = PROJ / "p3_qc_panel/output/qc_panel_239.json"
            if qc_panel_path.exists():
                qc_data = json.loads(qc_panel_path.read_text())
                gp = qc_data.get("panel_genes", qc_data.get("genes", []))
                if isinstance(gp[0], dict):
                    gp = [g["gene"] if isinstance(g, dict) else g for g in gp]
            else:
                gp = ["UXS1","PLOD1","MALAT1","C1D","KIF1B","COL1A1","FN1","VIM","ACTA2","MYH9",
                      "TPM1","TAGLN","CNN1","MYL9","ACTG2","DES","MYH11","LMOD1","PDLIM3","MYLK",
                      "ITGA8","SYNPO2","MRVI1","PTGIS","GUCY1A1","NEXN","PPP1R12B","SORBS1","FLNC","SPARCL1"]

            overlap = [g for g in gp if g in adata.var_names]
            print(f"\n  30-gene panel overlap: {len(overlap)}/{len(gp)} genes")
            if overlap:
                print(f"    Found: {overlap}")

            # Save
            results = {
                "dataset": "GSE240556",
                "n_cells_total": int(adata.n_obs),
                "n_genes_total": int(adata.n_vars),
                "median_genes_per_cell": float(np.median(adata.obs['n_genes'])),
                "median_counts_per_cell": float(np.median(adata.obs['n_counts'])),
                "panel_gene_overlap": overlap,
                "n_panel_genes_found": len(overlap),
                "qc_stats": {
                    "n_genes_median": float(np.median(adata.obs['n_genes'])),
                    "n_counts_median": float(np.median(adata.obs['n_counts'])),
                }
            }
            json.dump(results, open(OUT / "snrna_seq_analysis.json", "w"), indent=2)
            print("\nSaved to snrna_seq_analysis.json")
        else:
            print("  Could not find required 10x files (matrix, barcodes, features)")
    except ImportError as e:
        print(f"  Missing package: {e}")
        print("  Install: pip install anndata scanpy scipy")
    except Exception as e:
        print(f"  Processing error: {e}")
        import traceback
        traceback.print_exc()
else:
    print("\nNo extractable data found in supplementary files")
    print("These may be FTP directory listings, not data archives")
    print("Manual download from GEO may be needed")

print("=" * 60)
print("DONE")
