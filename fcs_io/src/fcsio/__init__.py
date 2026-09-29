"""fcsio: a stdlib reader for FCS 3.0/3.1 flow-cytometry files.

Reads the HEADER/TEXT/DATA layout directly — no pandas or numpy needed
for the parse itself. `read` returns (events, meta) where events is a
list of per-channel float lists and meta is the TEXT keyword dict.

Vendor quirks handled and tested:

  - TEXT delimiter is the first byte of the segment, not assumed to be
    the documented \\x0c — a doubled delimiter inside a value is an
    escaped literal.
  - Header DATA offsets may be zero; $BEGINDATA/$ENDDATA keywords are
    the fallback per spec.
  - $BYTEORD arbitrary byte orders for 4-byte values ("1,2,3,4" and
    "4,3,2,1" both seen; others supported when uniform).
  - $DATATYPE I/F/D; INTEGER requires uniform $PnB (mixed widths are a
    clear error, not silent misparsing).
  - $SPILLOVER extracted to a matrix keyed by channel order.
  - Multi-dataset files: only the first dataset is read; $NEXTDATA is
    surfaced in meta rather than silently dropping later datasets.

Deliberately absent: log amplification ($PnE) is metadata here — the
caller decides when to transform, matching cytof_qc's approach.
"""

from .reader import read, FCSError  # noqa: F401
from .spillover import spillover  # noqa: F401
from .writer import write  # noqa: F401
