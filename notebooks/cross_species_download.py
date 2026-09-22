"""Download cross-species datasets with rate-limit handling."""
import json
import time
import urllib.error
import urllib.request
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
PROJ = Path(__file__).resolve().parents[1]
CROSS = PROJ / "cross_species_validation"
CROSS.mkdir(exist_ok=True)

print("=" * 60)
print("CROSS-SPECIES DOWNLOAD (rate-limit aware)")
print("=" * 60)

def ncbi_request(url, retries=3, delay=3):
    """Make NCBI request with retries and rate-limit backoff."""
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "RaoLab/1.0 (barlowa124@gmail.com)"})
            with urllib.request.urlopen(req, timeout=60) as resp:
                return resp.read().decode("utf-8", errors="ignore")
        except urllib.error.HTTPError as e:
            if e.code == 429:
                wait = delay * (2 ** attempt)
                print(f"     Rate limited, waiting {wait}s...")
                time.sleep(wait)
            else:
                print(f"     HTTP {e.code}: {e}")
                return None
        except Exception as e:
            print(f"     Error: {e}")
            time.sleep(delay)
    return None

# ── 1. Fetch metadata for all datasets ────────────────────
print("\n1. Fetching dataset metadata...")

datasets = {
    "GSE173199": "Bovine satellite cell diff (cultivated meat)",
    "GSE240556": "Bovine snRNA-seq diff (Nature Comms 2024)",
    "GSE206914": "Porcine scRNA-seq myogenesis",
}

for gse, desc in datasets.items():
    print(f"\n   {gse} ({desc})...")
    time.sleep(1.5)  # Be nice to NCBI
    
    # Search for GEO ID
    search_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term={gse}&retmode=xml"
    text = ncbi_request(search_url)
    if not text:
        continue
    
    import xml.etree.ElementTree as ET
    root = ET.fromstring(text)
    ids = [e.text for e in root.findall(".//Id")]
    if not ids:
        print("     No results")
        continue
    
    geo_id = ids[0]
    time.sleep(1)
    
    # Get summary
    summary_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gds&id={geo_id}&retmode=xml"
    text = ncbi_request(summary_url)
    if not text:
        continue
    
    root = ET.fromstring(text)
    doc = root.find(".//DocSum")
    if doc is None:
        continue
    
    def get_item(name):
        el = doc.find(f".//Item[@Name='{name}']")
        return el.text if el is not None else "N/A"
    
    title = get_item("title")
    organism = get_item("taxon")
    n_samples = get_item("n_samples")
    summary = get_item("summary")
    ptype = get_item("ptype")
    
    print(f"     Title: {title[:120]}")
    print(f"     Organism: {organism}, Samples: {n_samples}, Type: {ptype}")
    
    # Get sample accessions
    time.sleep(1)
    samples_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=gds&db=sra&id={geo_id}&retmode=xml"
    text = ncbi_request(samples_url)
    sra_ids = []
    if text:
        root = ET.fromstring(text)
        sra_ids = [e.text for e in root.findall(".//Link/Id")]
    
    print(f"     SRA accessions: {len(sra_ids)}")
    
    # Save metadata
    meta = {
        "gse": gse, "geo_id": geo_id, "title": title,
        "organism": organism, "n_samples": n_samples, "type": ptype,
        "summary": (summary or "")[:500], "sra_ids": sra_ids[:10],
        "description": desc,
    }
    json.dump(meta, open(CROSS / f"{gse}_metadata.json", "w"), indent=2)

# ── 2. Download supplementary files ───────────────────────
print("\n\n2. Downloading supplementary files...")

for gse in datasets:
    print(f"\n   {gse}...")
    time.sleep(2)
    
    # Get the series page to find supplementary links
    url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse}"
    text = ncbi_request(url)
    if not text:
        continue
    
    import re
    # Find all download links
    suppl_urls = re.findall(r'https?://[^\s"\']+\.(?:txt|csv|gz|zip|tar)[^\s"\']*', text)
    # Also find FTP links
    suppl_urls += re.findall(r'ftp://[^\s"\']+', text)
    
    # Deduplicate
    suppl_urls = list(set(suppl_urls))
    
    out_dir = CROSS / gse
    out_dir.mkdir(exist_ok=True)
    
    for url in suppl_urls[:8]:
        fname = url.split("/")[-1].split("?")[0]
        if not fname:
            fname = f"suppl_{hash(url) % 10000}.txt"
        fpath = out_dir / fname
        
        if fpath.exists():
            print(f"     Already have: {fname}")
            continue
        
        try:
            time.sleep(1)
            req = urllib.request.Request(url, headers={"User-Agent": "RaoLab/1.0"})
            with urllib.request.urlopen(req, timeout=120) as resp:
                data = resp.read()
            fpath.write_bytes(data)
            print(f"     Downloaded: {fname} ({len(data)} bytes)")
        except Exception as e:
            print(f"     Failed: {fname} - {e}")

# ── 3. Summary ────────────────────────────────────────────
print("\n\n3. Files downloaded:")
for gse in datasets:
    d = CROSS / gse
    if d.exists():
        files = list(d.glob("*"))
        print(f"   {gse}: {len(files)} files")
        for f in files:
            print(f"     {f.name} ({f.stat().st_size} bytes)")

print(f"\nOutput: {CROSS}")
print("=" * 60)
print("DOWNLOAD COMPLETE")
print("=" * 60)