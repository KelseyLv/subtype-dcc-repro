"""
Batch-run repro/evaluate.py metrics for all nine TCGA cancers used in Subtype-DCC Table 1.

Usage (from repository root):

    .venv\\Scripts\\python.exe repro\\batch_evaluate.py

Optional:

    .venv\\Scripts\\python.exe repro\\batch_evaluate.py --strict
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

REPRO_DIR = Path(__file__).resolve().parent
REPO_ROOT = REPRO_DIR.parent
if str(REPRO_DIR) not in sys.path:
    sys.path.insert(0, str(REPRO_DIR))

from data_io import load_clinical_table  # noqa: E402
from evaluate import (  # noqa: E402
    clinical_significant_count,
    load_labels_dcc,
    merge_labels_clinical,
    survival_neg_log10_p,
)

# Briefings in Bioinformatics Table 1 cancers (Subtype-DCC paper)
DEFAULT_CANCERS: List[str] = [
    "BRCA",
    "BLCA",
    "KIRC",
    "LUAD",
    "PAAD",
    "SKCM",
    "STAD",
    "UCEC",
    "UVM",
]


def evaluate_labels_to_dict(
    cancer: str,
    labels_path: Path,
    clinical_path: Path,
) -> Dict[str, Any]:
    clin = load_clinical_table(clinical_path, cancer)
    labels = load_labels_dcc(labels_path)
    merged = merge_labels_clinical(labels, clin)
    if merged.empty:
        raise ValueError("No overlap between labels and clinical barcodes.")

    nlp = survival_neg_log10_p(merged)
    n_sig, pvals = clinical_significant_count(merged)
    return {
        "cancer": cancer,
        "method": None,
        "labels_path": str(labels_path.resolve()),
        "n_samples_eval": int(len(merged)),
        "neg_log10_logrank_p": nlp,
        "n_significant_clinical": n_sig,
        "clinical_pvals": pvals,
    }


def run_one(
    cancer: str,
    labels_path: Path,
    clinical_path: Path,
) -> Dict[str, Any]:
    """Backward-compatible alias (Subtype-DCC .dcc batch)."""
    return evaluate_labels_to_dict(cancer, labels_path, clinical_path)


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Batch evaluate Subtype-DCC *.dcc vs clinical (paper-style metrics)."
    )
    ap.add_argument(
        "--clinical",
        type=Path,
        default=REPO_ROOT / "data" / "clinical" / "cBioportal_data.tsv",
        help="Clinical TSV (GerkeLab or legacy PANCAN).",
    )
    ap.add_argument(
        "--labels-dir",
        type=Path,
        default=REPO_ROOT
        / "vendor"
        / "Subtype-DCC"
        / "Subtype-DCC"
        / "results",
        help="Directory containing <CANCER>.dcc files.",
    )
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=REPO_ROOT / "repro_out",
        help="Directory for per-cancer JSON and batch summary.",
    )
    ap.add_argument(
        "--cancers",
        type=str,
        default=",".join(DEFAULT_CANCERS),
        help="Comma-separated cancer codes (default: nine paper cancers).",
    )
    ap.add_argument(
        "--strict",
        action="store_true",
        help="Exit with code 1 if any cancer is missing .dcc or fails.",
    )
    args = ap.parse_args()

    cancers = [c.strip().upper() for c in args.cancers.split(",") if c.strip()]
    args.out_dir.mkdir(parents=True, exist_ok=True)

    summary: Dict[str, Any] = {
        "clinical_path": str(args.clinical.resolve()),
        "labels_dir": str(args.labels_dir.resolve()),
        "results": {},
    }

    failures = 0
    for cancer in cancers:
        labels_path = args.labels_dir / f"{cancer}.dcc"
        if not labels_path.is_file():
            summary["results"][cancer] = {
                "status": "missing_dcc",
                "labels_path": str(labels_path),
            }
            failures += 1
            print(f"[skip] {cancer}: missing {labels_path}", flush=True)
            continue
        try:
            out = run_one(cancer, labels_path, args.clinical)
            out["status"] = "ok"
            per_path = args.out_dir / f"{cancer}.metrics.json"
            per_path.write_text(
                json.dumps(out, indent=2, ensure_ascii=False),
                encoding="utf-8",
            )
            summary["results"][cancer] = {
                "status": "ok",
                "metrics_json": str(per_path.resolve()),
                "n_samples_eval": out["n_samples_eval"],
                "neg_log10_logrank_p": out["neg_log10_logrank_p"],
                "n_significant_clinical": out["n_significant_clinical"],
            }
            print(
                f"[ok] {cancer}: n={out['n_samples_eval']}, "
                f"-log10 logrank p={out['neg_log10_logrank_p']:.4f}, "
                f"clinical_sig={out['n_significant_clinical']}/6 -> {per_path}",
                flush=True,
            )
        except Exception as e:
            failures += 1
            summary["results"][cancer] = {"status": "error", "error": str(e)}
            print(f"[error] {cancer}: {e}", flush=True)

    summary_path = args.out_dir / "batch_metrics_summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    print(f"\nWrote summary: {summary_path}", flush=True)

    if args.strict and failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
