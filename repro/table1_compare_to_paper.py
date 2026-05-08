"""
Merge repro_out/table1_reproduction_metrics.json with data/paper_table1_reference.json.

Produces:
  repro_out/table1_reproduction_vs_paper.json  — flat rows + paper columns + deltas (if paper present)
  repro_out/table1_reproduction_vs_paper.md    — Markdown tables per method

Paper SNF row is bundled in data/paper_table1_reference.json; add other methods from the PDF when needed.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

REPO_ROOT = Path(__file__).resolve().parents[1]


def _paper_pair(
    paper: Dict[str, Any], method: str, cancer: str
) -> Optional[Tuple[float, int]]:
    block = paper.get(method)
    if not isinstance(block, dict):
        return None
    raw = block.get(cancer)
    if raw is None or not isinstance(raw, (list, tuple)) or len(raw) < 2:
        return None
    try:
        return float(raw[0]), int(raw[1])
    except (TypeError, ValueError):
        return None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--repro",
        type=Path,
        default=REPO_ROOT / "repro_out" / "table1_reproduction_metrics.json",
    )
    ap.add_argument(
        "--paper",
        type=Path,
        default=REPO_ROOT / "data" / "paper_table1_reference.json",
    )
    ap.add_argument(
        "--out-json",
        type=Path,
        default=REPO_ROOT / "repro_out" / "table1_reproduction_vs_paper.json",
    )
    ap.add_argument(
        "--out-md",
        type=Path,
        default=REPO_ROOT / "repro_out" / "table1_reproduction_vs_paper.md",
    )
    args = ap.parse_args()

    repro_data = json.loads(args.repro.read_text(encoding="utf-8"))
    paper_data = json.loads(args.paper.read_text(encoding="utf-8"))

    merged_rows: List[Dict[str, Any]] = []
    for row in repro_data.get("rows", []):
        if row.get("status") != "ok":
            merged_rows.append(dict(row))
            continue
        cancer = row["cancer"]
        method = row["method"]
        pair = _paper_pair(paper_data, method, cancer)
        base = {
            "cancer": cancer,
            "method": method,
            "n_samples_eval": row.get("n_samples_eval"),
            "repro_neg_log10_logrank_p": row.get("neg_log10_logrank_p"),
            "repro_n_significant_clinical": row.get("n_significant_clinical"),
            "labels_path": row.get("labels_path"),
            "status": "ok",
        }
        if pair is None:
            base["paper_neg_log10_logrank_p"] = None
            base["paper_n_significant_clinical"] = None
            base["delta_neg_log10_logrank_p"] = None
            base["delta_n_significant_clinical"] = None
            base["paper_available"] = False
        else:
            plp, pn = pair
            rlp = float(row["neg_log10_logrank_p"])
            rn = int(row["n_significant_clinical"])
            base["paper_neg_log10_logrank_p"] = plp
            base["paper_n_significant_clinical"] = pn
            base["delta_neg_log10_logrank_p"] = round(rlp - plp, 6)
            base["delta_n_significant_clinical"] = rn - pn
            base["paper_available"] = True
        merged_rows.append(base)

    out_obj = {
        "citation": paper_data.get("citation"),
        "note": paper_data.get("note"),
        "repro_source": str(args.repro.resolve()),
        "paper_source": str(args.paper.resolve()),
        "rows": merged_rows,
    }
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(out_obj, indent=2, ensure_ascii=False), encoding="utf-8")

    # Markdown: one table per method
    by_method: Dict[str, List[Dict[str, Any]]] = {}
    for r in merged_rows:
        if r.get("status") != "ok":
            continue
        by_method.setdefault(r["method"], []).append(r)

    lines = [
        "# Table 1 reproduction vs paper",
        "",
        f"Reproduction metrics from `{args.repro.as_posix()}`; paper reference from `{args.paper.as_posix()}`.",
        "",
        "Columns: **repro** = this repository (`cBioportal_data.tsv` + `evaluate.py`); **paper** = bundled Table 1 values where present; **Δ** = repro − paper.",
        "",
    ]
    method_order = ["subtype_dcc", "snf", "kmeans", "spectral", "nmf"]
    ordered = [m for m in method_order if m in by_method] + sorted(
        m for m in by_method if m not in method_order
    )
    for method in ordered:
        lines.append(f"## {method}")
        lines.append("")
        lines.append(
            "| Cancer | n (repro) | repro −log10 p | repro clin | paper −log10 p | paper clin | Δ −log10 p | Δ clin |"
        )
        lines.append("|--------|-----------|----------------|-------------|-----------------|------------|------------|--------|")
        for r in sorted(by_method[method], key=lambda x: x["cancer"]):
            pp = r.get("paper_neg_log10_logrank_p")
            pc = r.get("paper_n_significant_clinical")
            ps = f"{pp:.2f}" if pp is not None and not (isinstance(pp, float) and math.isnan(pp)) else "—"
            cs = str(pc) if pc is not None else "—"
            dlp = r.get("delta_neg_log10_logrank_p")
            dlc = r.get("delta_n_significant_clinical")
            ds = f"{dlp:+.2f}" if dlp is not None else "—"
            dcs = f"{dlc:+d}" if dlc is not None else "—"
            lines.append(
                f"| {r['cancer']} | {r['n_samples_eval']} | "
                f"{r['repro_neg_log10_logrank_p']:.4f} | {r['repro_n_significant_clinical']} | "
                f"{ps} | {cs} | {ds} | {dcs} |"
            )
        lines.append("")

    args.out_md.write_text("\n".join(lines), encoding="utf-8")
    print(args.out_json)
    print(args.out_md)


if __name__ == "__main__":
    main()
