"""Download and process porcine dataset GSE206914."""
import re
import time
import urllib.request
from pathlib import Path

CROSS = Path(r"C:\Users\asdf\CascadeProjects\rao_lab_ml\cultivated_meat_projects\cross_species_validation")
OUT = CROSS / "GSE206914"
OUT.mkdir(exist_ok=True)

print("=" * 60)
print("PORCINE DATASET GSE206914")
print("=" * 60)

# Check sub-series
for gse in ['GSE206914', 'GSE206912', 'GSE206913']:
    print(f"\n{gse}:")
    try:
        req = urllib.request.Request(
            f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse}",
            headers={"User-Agent": "RaoLab/1.0"}
        )
        with urllib.request.urlopen(req, timeout=30) as r:
            text = r.read().decode("utf-8", errors="ignore")
        
        # Title
        title_m = re.search(r'<td[^>]*>Title</td>.*?<td[^>]*>(.*?)</td>', text, re.DOTALL)
        if title_m:
            t = re.sub(r'<[^>]+>', '', title_m.group(1)).strip()
            print(f"  Title: {t[:150]}")
        
        # Organism
        org_m = re.search(r'<td[^>]*>Organism</td>.*?<td[^>]*>(.*?)</td>', text, re.DOTALL)
        if org_m:
            print(f"  Organism: {re.sub(r'<[^>]+>', '', org_m.group(1)).strip()}")
        
        # Samples
        samp_m = re.findall(r'GSM\d+', text)
        print(f"  GSM samples: {len(set(samp_m))}")
        
        # Supplementary
        suppl = re.findall(r'https?://[^\s\"\'<>]+\.(?:txt|gz|csv|tsv|zip|mtx|h5|h5ad|rds)[^\s\"\'<>]*', text)
        for u in suppl[:5]:
            print(f"  Suppl: {u[-100:]}")
        
        # Matrix
        matrix = re.findall(r'/geo/download/[^\s\"\'<>]+', text)
        for m in matrix[:3]:
            print(f"  Matrix: https://www.ncbi.nlm.nih.gov{m}")
            
    except Exception as e:
        print(f"  Error: {e}")
    time.sleep(1.5)

print("\nDone checking porcine datasets.")
print("=" * 60)