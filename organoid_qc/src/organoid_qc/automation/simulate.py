"""Validate the Flex protocol in software via opentrons.simulate.

`run_simulation` executes the protocol exactly as the robot's software
would (same Protocol API engine), returning the full run log, every
aspirate/dispense is an inspectable step, so volume accounting can be
asserted in tests.
"""

from __future__ import annotations

import io
import json
from pathlib import Path

PROTOCOL = Path(__file__).with_name("flex_organoid_dosing.py")


def run_simulation() -> list[dict]:
    """Simulate the dosing protocol; return the run-log entries."""
    from opentrons.simulate import simulate

    with open(PROTOCOL, "rb") as f:
        runlog, _bundle = simulate(
            f, custom_labware_paths=[], propagate_logs=False
        )
    return runlog


def runlog_text(entries: list[dict]) -> str:
    buf = io.StringIO()
    for e in entries:
        payload = e.get("payload", {})
        buf.write(
            f"{e.get('level', 'info')}: {payload.get('text', e)}\n"
        )
    return buf.getvalue()


if __name__ == "__main__":
    entries = run_simulation()
    print(json.dumps(entries, indent=2, default=str)[:4000])
    print(f"\n{len(entries)} logged steps")
