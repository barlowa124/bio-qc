"""Download GSE240556 bovine snRNA-seq (Nature Comms 2024, cultivated meat heterogeneity)."""
import gzip
import re
import shutil
import time
import urllib.request
from pathlib import Path

CROSS = Path(r"C:\Users\asdf\CascadeProjects\rao_lab_ml\cultivated_meat_projects\cross_species_validation")
OUT = CROSS / "GSE240556"
OUT.mkdir(exist_ok=True)

print("=" * 60)
print("DOWNLOADING GSE240556 — Bovine snRNA-seq (Nature Comms 2024)")
print("=" * 60)

# 1. Get the series page
print("\n1. Fetching series page...")
url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE240556"
req = urllib.request.Request(url, headers={"User-Agent": "RaoLab/1.0"})
with urllib.request.urlopen(req, timeout=30) as r:
    text = r.read().decode("utf-8", errors="ignore")

# 2. Find supplementary files
print("\n2. Finding supplementary files...")
suppl_urls = re.findall(r'https?://[^\s\"\'<>]+\.(?:txt|gz|csv|tsv|zip|mtx|h5|h5ad|rds|tar)[^\s\"\'<>]*', text)
matrix_urls = re.findall(r'/geo/download/[^\s\"\'<>]+', text)

print(f"   Supplementary URLs: {len(suppl_urls)}")
for u in suppl_urls[:10]:
    print(f"     {u[-100:]}")

print(f"   Matrix links: {len(matrix_urls)}")
for m in matrix_urls[:5]:
    full = f"https://www.ncbi.nlm.nih.gov{m}"
    print(f"     {full}")

# 3. Download supplementary files
print("\n3. Downloading files...")
all_urls = suppl_urls + [f"https://www.ncbi.nlm.nih.gov{m}" for m in matrix_urls]

for url in all_urls[:10]:
    fname = url.split("/")[-1].split("?")[0]
    if not fname or len(fname) < 3:
        fname = f"suppl_{abs(hash(url)) % 100000}.txt"
    fpath = OUT / fname
    
    if fpath.exists():
        print(f"   Already have: {fname} ({fpath.stat().st_size} bytes)")
        continue
    
    try:
        time.sleep(1)
        req = urllib.request.Request(url, headers={"User-Agent": "RaoLab/1.0"})
        with urllib.request.urlopen(req, timeout=120) as r:
            data = r.read()
        fpath.write_bytes(data)
        print(f"   Downloaded: {fname} ({len(data)} bytes)")
        
        # Decompress if gzipped
        if fname.endswith('.gz'):
            csv_path = OUT / fname.replace('.gz', '')
            try:
                with gzip.open(fpath, 'rb') as f_in:
                    with open(csv_path, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                print(f"     Decompressed: {csv_path.stat().st_size} bytes")
            except (OSError, gzip.BadGzipFile):
                pass
    except Exception as e:
        print(f"   Failed: {fname} - {e}")

# 4. Summary
print(f"\n4. Files in {OUT}:")
for f in sorted(OUT.glob("*")):
    print(f"   {f.name} ({f.stat().st_size} bytes)")

print("=" * 60)
print("GSE240556 DOWNLOAD COMPLETE")
print("=" * 60)