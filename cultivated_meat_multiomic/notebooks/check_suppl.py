import re
import time
import urllib.request
from pathlib import Path

CROSS = Path(r"C:\Users\asdf\CascadeProjects\rao_lab_ml\cultivated_meat_projects\cross_species_validation")

# Check sub-series too
targets = {
    "GSE173196": "Bovine SF_Diff subseries",
    "GSE173198": "Bovine Timecourse subseries", 
    "GSE240556": "Bovine snRNA-seq",
    "GSE206914": "Porcine scRNA-seq",
}

for gse, desc in targets.items():
    print(f"\n{gse} ({desc}):")
    try:
        req = urllib.request.Request(
            f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse}",
            headers={"User-Agent": "RaoLab/1.0"}
        )
        with urllib.request.urlopen(req, timeout=30) as r:
            text = r.read().decode("utf-8", errors="ignore")
        
        # Find supplementary file links
        urls = re.findall(r'https?://[^\s\"\'<>]+\.(?:txt|gz|csv|tsv|zip|xlsx)[^\s\"\'<>]*', text)
        for u in urls[:8]:
            print(f"  {u}")
        if not urls:
            print("  No supplementary files found on page")
            
        # Also check for "Series Matrix File(s)" link
        matrix = re.findall(r'/geo/download/[^\s\"\'<>]+', text)
        for m in matrix[:3]:
            print(f"  Matrix: https://www.ncbi.nlm.nih.gov{m}")
            
    except Exception as e:
        print(f"  Error: {e}")
    time.sleep(1.5)

print("\nDone.")