"""Synthetic FCS generator for the QC pipeline demo.

Writes small gated-population FCS files (Levine-style
``<sample>_<population>_cells.fcs``) via fcsio, so the pipeline runs end
to end without downloading the Levine benchmark. Channels match the
cytof_qc conventions: phenotypic markers plus QC-only channels
(Time/DNA/Viability) that must be excluded from clustering.

Deterministic seed — demo files are fixtures, not biological findings.
"""

from __future__ import annotations

import random
from pathlib import Path

CHANNELS = ["Time", "DNA1", "DNA2", "Viability",
            "CD45", "CD3", "CD4", "CD8", "CD19", "CD56", "CD11c", "CD14"]

# Population signatures: (marker, mean) pairs on arcsinh-scale-ish values.
POPULATIONS = {
    "CD4 T": {"CD45": 5.0, "CD3": 4.5, "CD4": 4.2},
    "CD8 T": {"CD45": 5.0, "CD3": 4.5, "CD8": 4.0},
    "B": {"CD45": 5.0, "CD19": 4.8},
    "NK": {"CD45": 5.0, "CD56": 4.4},
    "Monocytes": {"CD45": 5.0, "CD11c": 4.0, "CD14": 4.6},
    "NotGated": {},                       # residual population, no signal
}


def main(out_dir: str, n_events: int = 400, seed: int = 0) -> None:
    from fcsio import write
    rng = random.Random(seed)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    for pop, means in POPULATIONS.items():
        events = []
        for i in range(n_events):
            row = []
            for ch in CHANNELS:
                if ch == "Time":
                    row.append(float(i))
                elif ch.startswith(("DNA", "Viability")):
                    row.append(rng.gauss(5.0, 0.3))
                elif ch in means:
                    row.append(rng.gauss(means[ch], 0.5))
                else:
                    row.append(rng.gauss(0.5, 0.4))     # background
            events.append(row)
        write(str(out / f"Demo_{pop}_cells.fcs"),
              events, {"channels": CHANNELS})


if __name__ == "__main__":
    import sys
    main(sys.argv[1] if len(sys.argv) > 1 else "data/demo",
         int(sys.argv[2]) if len(sys.argv) > 2 else 400)
