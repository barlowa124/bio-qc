"""Verified-claims layer: every number in the markdown report must bind
to a recorded value in the run's JSON artifacts.

Ported from the oncology_coscientist draft verifier (trust-tools). The
pattern: flatten every numeric leaf of results/*.json into a value pool,
extract every numeric token from the report text, and require each token
to sit within display-rounding tolerance of a pooled value. A claim that
no recorded computation produced is flagged as unbound.

Cluster names are exempt like metric-key integers are in oncocs: "3" in
"cluster 3" is a label, not a measurement.
"""

from __future__ import annotations

import re
from typing import Any

FORBIDDEN_PHRASES = [
    "clinically validated",
    "proves",
    "causes",
    "state-of-the-art",
]

# Same boundary logic as oncocs: a '-' only attaches when not preceded by
# a word char/digit/dot, so "0.5-0.6" yields two tokens while " -0.5"
# keeps its sign.
_NUM_RE = re.compile(r"(?<![\w.%])-?\d(?:[\d,]*\d)?(?:\.\d+)?\s*%?(?![\w.%])")


def flatten_results(obj: Any, prefix: str = "") -> dict[str, float]:
    """Flatten every numeric leaf to {dotted.path: value}."""
    out: dict[str, float] = {}
    if isinstance(obj, dict):
        for k, v in obj.items():
            out.update(flatten_results(v, f"{prefix}{k}."))
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            out.update(flatten_results(v, f"{prefix}{i}."))
    elif isinstance(obj, bool):
        pass
    elif isinstance(obj, (int, float)):
        out[prefix[:-1]] = float(obj)
    return out


def extract_numbers(text: str) -> list[dict]:
    """Numeric tokens with display precision and context."""
    tokens = []
    for m in _NUM_RE.finditer(text):
        raw = m.group(0).strip()
        is_pct = raw.endswith("%")
        val = raw.rstrip("%").replace(",", "")
        decimals = len(val.split(".")[1]) if "." in val else 0
        ctx = text[max(0, m.start() - 30): m.end() + 30].replace("\n", " ")
        tokens.append({"token": raw, "value": float(val), "is_pct": is_pct,
                       "decimals": decimals, "context": ctx.strip(),
                       "pos": m.start()})
    return tokens


def verify_markdown(text: str, flat_values: dict[str, float],
                    labels: set[str] | None = None) -> dict:
    """Bind every numeric token in `text` to a pooled recorded value.

    `labels` holds exempt identifiers (cluster names): a zero-decimal
    token equal to a label is an identifier, not a measurement.
    """
    values = list(flat_values.values())
    labels = {str(x) for x in (labels or set())}
    unbound, forbidden = [], []

    for tok in extract_numbers(text):
        v = tok["value"]
        tol = 0.5 * 10 ** (-tok["decimals"])
        # a % token also matches the fraction form of the same value
        ok = any(abs(v - fv) <= tol for fv in values)
        if not ok and tok["is_pct"]:
            ok = any(abs(v / 100.0 - fv) <= tol / 100.0 for fv in values)
        if not ok and tok["decimals"] == 0 and str(int(v)) in labels:
            ok = True
        if not ok:
            unbound.append({"token": tok["token"], "context": tok["context"]})

    low = text.lower()
    for phrase in FORBIDDEN_PHRASES:
        if re.search(r"\b" + re.escape(phrase) + r"\b", low):
            forbidden.append(phrase)

    return {"passed": not unbound and not forbidden,
            "n_claims": len(extract_numbers(text)),
            "unbound_claims": unbound,
            "forbidden": forbidden}
