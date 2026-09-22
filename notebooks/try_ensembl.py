"""Try Ensembl REST API with small batches."""
import json
import urllib.request
from pathlib import Path

CROSS = Path(r"C:\Users\asdf\CascadeProjects\rao_lab_ml\cultivated_meat_projects\cross_species_validation")

# Load bovine gene IDs
import pandas as pd

sf = pd.read_csv(CROSS / "GSE173199/sf_diff_counts.csv", index_col=0, nrows=100)
ids = [g.split(".")[0] for g in sf.index.tolist()]
print(f"Testing with {len(ids)} IDs")

# Try POST with small batch
batch = ids[:10]
payload = json.dumps({"ids": batch})
print(f"Payload: {payload[:100]}...")

try:
    req = urllib.request.Request(
        "https://rest.ensembl.org/lookup/id",
        data=payload.encode(),
        headers={"Content-Type": "application/json", "Accept": "application/json"}
    )
    with urllib.request.urlopen(req, timeout=30) as r:
        data = json.loads(r.read())
        print(f"Got {len(data)} results")
        for eid, info in list(data.items())[:5]:
            print(f"  {eid} -> {info.get('display_name', 'N/A')}")
except Exception as e:
    print(f"Error: {type(e).__name__}: {e}")
    # Try alternative: requests library
    try:
        import requests
        print("\nTrying with requests library...")
        r = requests.post(
            "https://rest.ensembl.org/lookup/id",
            headers={"Content-Type": "application/json", "Accept": "application/json"},
            data=payload,
            timeout=30
        )
        print(f"Status: {r.status_code}")
        data = r.json()
        for eid, info in list(data.items())[:5]:
            print(f"  {eid} -> {info.get('display_name', 'N/A')}")
    except ImportError:
        print("requests not installed")
    except Exception as e2:
        print(f"requests error: {e2}")