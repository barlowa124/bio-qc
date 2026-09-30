"""Map all bovine ENSBTAG IDs to gene symbols via Ensembl REST."""
import json
import time
from pathlib import Path

import pandas as pd
import requests

CROSS = Path(r"C:\Users\asdf\CascadeProjects\rao_lab_ml\cultivated_meat_projects\cross_species_validation")
MAPPING_FILE = CROSS / "bovine_ensembl_to_symbol.json"

print("=" * 60)
print("MAPPING BOVINE ENSEMBL -> GENE SYMBOLS")
print("=" * 60)

# Load all bovine gene IDs
sf = pd.read_csv(CROSS / "GSE173199/sf_diff_counts.csv", index_col=0)
tc = pd.read_csv(CROSS / "GSE173199/timecourse_counts.csv", index_col=0)
all_ids = list(set([g.split(".")[0] for g in list(sf.index) + list(tc.index)]))
print(f"\nUnique Ensembl IDs: {len(all_ids)}")

# Check if already mapped
if MAPPING_FILE.exists():
    existing = json.loads(MAPPING_FILE.read_text())
    print(f"Existing mappings: {len(existing)}")
    remaining = [eid for eid in all_ids if eid not in existing]
    symbol_map = existing
else:
    remaining = all_ids
    symbol_map = {}

print(f"Remaining to map: {len(remaining)}")

# Batch query
BATCH = 500
for i in range(0, len(remaining), BATCH):
    batch = remaining[i:i+BATCH]
    payload = json.dumps({"ids": batch})
    
    try:
        r = requests.post(
            "https://rest.ensembl.org/lookup/id",
            headers={"Content-Type": "application/json", "Accept": "application/json"},
            data=payload,
            timeout=60
        )
        if r.status_code == 200:
            data = r.json()
            for eid, info in data.items():
                if info is not None and isinstance(info, dict):
                    symbol_map[eid] = info.get("display_name", eid)
                else:
                    symbol_map[eid] = eid
            print(f"  Batch {i//BATCH + 1}: {len(batch)} sent, {len(symbol_map)} total mapped")
        elif r.status_code == 429:
            wait = 5
            print(f"  Rate limited, waiting {wait}s...")
            time.sleep(wait)
            # Retry this batch
            r = requests.post(
                "https://rest.ensembl.org/lookup/id",
                headers={"Content-Type": "application/json", "Accept": "application/json"},
                data=payload,
                timeout=60
            )
            if r.status_code == 200:
                data = r.json()
                for eid, info in data.items():
                    if info is not None and isinstance(info, dict):
                        symbol_map[eid] = info.get("display_name", eid)
                    else:
                        symbol_map[eid] = eid
        else:
            print(f"  Batch {i//BATCH + 1}: HTTP {r.status_code}")
    except Exception as e:
        print(f"  Batch {i//BATCH + 1} error: {e}")
    
    # Save progress
    MAPPING_FILE.write_text(json.dumps(symbol_map))
    time.sleep(0.3)

print(f"\nFinal mappings: {len(symbol_map)}")
# Show some examples
for eid in list(symbol_map.keys())[:10]:
    print(f"  {eid} -> {symbol_map[eid]}")

print(f"\nSaved to {MAPPING_FILE}")
print("=" * 60)
print("MAPPING COMPLETE")
print("=" * 60)