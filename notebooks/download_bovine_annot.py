"""Map ENSBTAG to symbols using local GTF download."""
import gzip
import time
import urllib.request
from pathlib import Path

CROSS = Path(r"C:\Users\asdf\CascadeProjects\rao_lab_ml\cultivated_meat_projects\cross_species_validation")

# Try multiple Ensembl FTP URLs
urls = [
    "https://ftp.ensembl.org/pub/release-113/tsv/bos_taurus/Bos_taurus.ARS-UCD2.0.113.gene.txt.gz",
    "https://ftp.ensembl.org/pub/release-112/tsv/bos_taurus/Bos_taurus.ARS-UCD2.0.112.gene.txt.gz",
    "https://ftp.ensembl.org/pub/release-111/tsv/bos_taurus/Bos_taurus.ARS-UCD2.0.111.gene.txt.gz",
    "https://ftp.ensembl.org/pub/release-110/tsv/bos_taurus/Bos_taurus.ARS-UCD2.0.110.gene.txt.gz",
]

out_path = CROSS / "bovine_genes.tsv.gz"

for url in urls:
    print(f"Trying: {url}")
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "RaoLab/1.0"})
        with urllib.request.urlopen(req, timeout=60) as r:
            data = r.read()
            out_path.write_bytes(data)
            print(f"  SUCCESS: {len(data)} bytes")
            
            # Parse it
            with gzip.open(out_path, 'rt') as f:
                header = f.readline().strip()
                print(f"  Header: {header[:100]}")
                for i, line in enumerate(f):
                    if i < 3:
                        parts = line.strip().split('\t')
                        print(f"  Row {i}: gene_id={parts[0]}, symbol={parts[2] if len(parts) > 2 else 'N/A'}")
                    else:
                        break
            break
    except Exception as e:
        print(f"  Failed: {e}")
    time.sleep(1)

print("\nDone.")