import re

import pytest

pytest.importorskip("opentrons")

from organoid_qc.automation.simulate import run_simulation, runlog_text


@pytest.fixture(scope="module")
def runlog():
    return runlog_text(run_simulation())


def _dispenses(runlog, labware):
    return [
        line
        for line in runlog.splitlines()
        if "Dispensing" in line and labware in line
    ]


def test_protocol_simulates_end_to_end(runlog):
    assert "Protocol complete" in runlog
    assert "Error" not in runlog


def test_dosing_transfers_count_and_volume(runlog):
    doses = _dispenses(runlog, "organoid plate")
    # 8 dilution points x 3 replicates = 24 doses of 10 uL
    assert len(doses) == 24
    assert all("10.0 uL" in d for d in doses)


def test_dilution_series_filled(runlog):
    fills = [
        line
        for line in runlog.splitlines()
        if "Transferring 200.0" in line
        and "reagent reservoir" in line
        and "serial dilution block" in line
    ]
    assert len(fills) == 8  # eight series wells pre-filled


def test_control_columns_untouched(runlog):
    doses = _dispenses(runlog, "organoid plate")
    wells = {
        re.search(r"into ([A-H]\d+) of organoid plate", d).group(1)
        for d in doses
    }
    # only columns 1-8 (one per dilution point) receive compound
    assert all(int(w[1:]) <= 8 for w in wells)
    assert len(wells) == 24
