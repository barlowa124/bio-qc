"""Opentrons Flex protocol: serial-dilution dosing of an organoid plate.

An 8-point, 3-fold compound dilution series is prepared in a deep-well
block and dosed in triplicate onto an organoid culture plate, leaving
columns 9-12 as vehicle-only controls. Volumes are chosen for a Flex
1-channel P1000; every transfer is logged by the Protocol API so the
simulation run can be audited step-by-step.

This is lab-automation code, not a biology claim: the same file runs
on the physical Flex, and `simulate.py` validates it in software.
"""

from opentrons import protocol_api

metadata = {
    "protocolName": "Organoid compound dosing (8pt 3-fold, triplicate)",
    "author": "organoid-qc",
    "description": (
        "Serial-dilute a test compound 3-fold over 8 points and dose a "
        "96-well organoid plate in triplicate with vehicle controls."
    ),
}
requirements = {"robotType": "Flex", "apiLevel": "2.19"}

DILUTION_FACTOR = 3.0
N_POINTS = 8
DILUENT_UL = 200.0
TRANSFER_UL = 100.0  # 100 into 200 -> 3-fold per step
DOSE_UL = 10.0
N_REPLICATES = 3
TOP_UP_MIX_REPS = 3


def run(ctx: protocol_api.ProtocolContext) -> None:
    dilution_block = ctx.load_labware(
        "nest_96_wellplate_2ml_deep", "C1", "serial dilution block"
    )
    organoid_plate = ctx.load_labware(
        "corning_96_wellplate_360ul_flat", "C2", "organoid plate"
    )
    reservoir = ctx.load_labware(
        "nest_12_reservoir_15ml", "C3", "reagent reservoir"
    )
    tiprack = ctx.load_labware("opentrons_flex_96_tiprack_1000ul", "D1")
    ctx.load_trash_bin("A3")
    p300 = ctx.load_instrument(
        "flex_1channel_1000", "right", tip_racks=[tiprack]
    )

    compound_stock = reservoir.wells_by_name()["A1"]
    diluent = reservoir.wells_by_name()["A2"]
    series = dilution_block.columns()[0][:N_POINTS]

    ctx.comment("Pre-fill dilution series with 200 uL media")
    p300.pick_up_tip()
    for well in series:
        p300.transfer(
            DILUENT_UL, diluent, well, new_tip="never", blow_out=True
        )
    p300.drop_tip()

    ctx.comment("Serial dilute compound 3-fold across 8 points")
    p300.transfer(
        TRANSFER_UL * DILUTION_FACTOR,
        compound_stock,
        series[0],
        mix_after=(TOP_UP_MIX_REPS, 200),
    )
    for src, dst in zip(series[:-1], series[1:]):
        p300.transfer(
            TRANSFER_UL, src, dst, mix_after=(TOP_UP_MIX_REPS, 200)
        )

    ctx.comment(
        f"Dose {DOSE_UL} uL in {N_REPLICATES} replicates (rows A-C); "
        "columns 9-12 left vehicle-only as controls"
    )
    for point, src in enumerate(series):
        for rep in range(N_REPLICATES):
            dest = organoid_plate.rows()[rep][point]
            p300.transfer(
                DOSE_UL, src, dest, new_tip="always", blow_out=True
            )
    ctx.comment("Protocol complete")
