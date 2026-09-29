# fcs_io

A dependency-free FCS 3.0/3.1 reader. `read(path)` returns `(events,
meta)`, a list of per-channel value rows plus the TEXT keyword dict with
derived `channels`, `stains`, `n_events`, `nextdata`.

```python
from fcsio import read, write, spillover
events, meta = read("sample.fcs")
names, matrix = spillover(meta)   # or None
write("out.fcs", events, {"channels": meta["channels"]})
```

Exists because cytof_qc leans on `fcsparser` for a binary format whose
quirks are worth owning. This reader handles the documented vendor
surface directly and fails loudly on what it does not support. It is
now wired in as the fallback reader in `cytof_qc/io.load_fcs`, so the
pipeline runs with no external FCS dependency. `write()` serializes
minimal valid FCS 3.0 files and exists for fixtures and demo data. The
`cytof_qc` workflow's `demo_data` rule uses it to generate the input
set, which the same package then reads back.

Quirks handled: delimiter-escaped values, header offsets blank with
$BEGINDATA/$ENDDATA fallback, little- and big-endian $BYTEORD,
$DATATYPE I/F/D, $SPILLOVER matrix extraction, $NEXTDATA surfaced.
Stated limits: INTEGER requires uniform $PnB (mixed widths
raise), only the first dataset of multi-set files is read, and no log
amplification is applied. $PnE stays metadata.

Tests build minimal FCS files in-memory (the fixture writer is part of
the evidence the format is understood). Real-file validation:
`scripts/crosscheck_fcsparser.py` compared every event value and
channel name against `fcsparser` on all 55 Levine_13dim/Levine_32dim
benchmark files — 432,671 events, bit-exact
(`validation/crosscheck_levine.json`). Rerun it after fetching the data
with `cytof_qc/scripts/fetch_data.sh`.
