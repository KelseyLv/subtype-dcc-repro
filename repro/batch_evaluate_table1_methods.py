"""
Evaluate Subtype-DCC + SNF + sklearn baselines (kmeans, spectral, nmf) for nine paper cancers.

Writes per-run JSON under repro_out/table1_eval/ and a flat summary
repro_out/table1_reproduction_metrics.json (45 rows when all label files exist).

Usage (repo root):

    .venv\\Scripts\\python.exe repro\\batch_evaluate_table1_methods.py
    .venv\\Scripts\\python.exe repro\\batch_evaluate_table1_methods.py --strict
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

REPRO_DIR = Path(__file__).resolve().parent
REPO_ROOT = REPRO_DIR.parent
if str(REPRO_DIR) not in sys.path:
    sys.path.insert(0, str(REPRO_DIR))

from batch_evaluate import DEFAULT_CANCERS, evaluate_labels_to_dict  # noqa: E402


METHOD_SPECS: List[Tuple[str, str]] = [
    ("subtype_dcc", "vendor/Subtype-DCC/Subtype-DCC/results/{cancer}.dcc"),
    ("snf", "repro_out/{cancer}.snf.tsv"),
    ("kmeans", "repro_out/{cancer}.kmeans.tsv"),
    ("spectral", "repro_out/{cancer}.spectral.tsv"),
    ("nmf", "repro_out/{cancer}.nmf.tsv"),
]


def labels_path(repo: Path, cancer: str, method: str) -> Path:
    rel = next(p for m, p in METHOD_SPECS if m == method)
    return repo / rel.format(cancer=cancer)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--clinical",
        type=Path,
        default=REPO_ROOT / "data" / "clinical" / "cBioportal_data.tsv",
    )
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=REPO_ROOT / "repro_out",
    )
    ap.add_argument(
        "--cancers",
        type=str,
        default=",".join(DEFAULT_CANCERS),
    )
    ap.add_argument("--strict", action="store_true")
    args = ap.parse_args()

    cancers = [c.strip().upper() for c in args.cancers.split(",") if c.strip()]
    eval_dir = args.out_dir / "table1_eval"
    eval_dir.mkdir(parents=True, exist_ok=True)

    rows: List[Dict[str, Any]] = []
    failures = 0

    for cancer in cancers:
        for method, _ in METHOD_SPECS:
            lp = labels_path(REPO_ROOT, cancer, method)
            if not lp.is_file():
                rows.append(
                    {
                        "cancer": cancer,
                        "method": method,
                        "status": "missing_labels",
                        "labels_path": str(lp),
                    }
                )
                failures += 1
                print(f"[skip] {cancer} {method}: missing {lp}", flush=True)
                continue
            try:
                out = evaluate_labels_to_dict(cancer, lp, args.clinical)
                out["method"] = method
                out["status"] = "ok"
                per = eval_dir / f"{cancer}_{method}.json"
                # Slim file without full pval dict for directory listing (optional full in summary)
                slim = {
                    "cancer": out["cancer"],
                    "method": method,
                    "labels_path": out["labels_path"],
                    "n_samples_eval": out["n_samples_eval"],
                    "neg_log10_logrank_p": out["neg_log10_logrank_p"],
                    "n_significant_clinical": out["n_significant_clinical"],
                    "status": "ok",
                }
                per.write_text(json.dumps(slim, indent=2, ensure_ascii=False), encoding="utf-8")
                rows.append(out)
                print(
                    f"[ok] {cancer} {method}: n={out['n_samples_eval']}, "
                    f"-log10p={out['neg_log10_logrank_p']:.4f}, clin={out['n_significant_clinical']}/6",
                    flush=True,
                )
            except Exception as e:
                failures += 1
                rows.append(
                    {
                        "cancer": cancer,
                        "method": method,
                        "status": "error",
                        "error": str(e),
                        "labels_path": str(lp),
                    }
                )
                print(f"[error] {cancer} {method}: {e}", flush=True)

    summary_path = args.out_dir / "table1_reproduction_metrics.json"
    summary_path.write_text(
        json.dumps(
            {
                "clinical_path": str(args.clinical.resolve()),
                "n_rows": len(rows),
                "rows": rows,
            },
            indent=2,
            ensure_ascii=False,
        ),
        encoding="utf-8",
    )
    print(f"\nWrote {summary_path}", flush=True)
    if args.strict and failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
