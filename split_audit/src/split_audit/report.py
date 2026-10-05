"""Render an audit result as a short markdown report."""
from __future__ import annotations


def render(name: str, result: dict, kind: str) -> str:
    lines = [f"# split-audit: {name}", ""]
    if kind == "scaffolds":
        lines += [
            f"- rows: {result['n_rows']}, scaffold ids: {result['n_scaffolds']}",
            f"- scaffold ids spanning the split boundary: "
            f"{result['n_shared_scaffolds']}",
            f"- leaked items: {result['n_leaked_items']}",
            f"- scaffold-disjoint: {result['scaffold_disjoint']}",
        ]
        if result["shared_scaffolds"]:
            lines.append(
                f"- shared scaffolds: {', '.join(result['shared_scaffolds'])}")
    else:
        lines += [
            f"- rows: {result['n_rows']}, splits: "
            f"{', '.join(result['splits'])}",
            f"- k-mer Jaccard threshold: {result['jaccard_threshold']} "
            f"(k={result['k']})",
            f"- cross-split pairs scored: {result['n_cross_split_pairs']}",
            f"- flagged pairs: {result['n_flagged_pairs']}",
            f"- flagged items: {result['n_flagged_items']}",
            f"- leakage-free: {result['leakage_free']}",
        ]
        if result["flagged_pairs"]:
            lines.append("- worst cross-split pairs:")
            for pair in result["flagged_pairs"]:
                lines.append(
                    f"  - {pair['a']} ~ {pair['b']}: {pair['similarity']:.3f}")
    return "\n".join(lines) + "\n"
