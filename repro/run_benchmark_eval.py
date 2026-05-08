"""
Batch-evaluate label files and write summary JSON for bar plots (Figure 2 style).
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from data_io import load_clinical_table
from evaluate import clinical_significant_count, merge_labels_clinical, survival_neg_log10_p


def load_any_labels(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    if "dcc" in df.columns and "cluster" not in df.columns:
        df = df.rename(columns={"dcc": "cluster"})
    if "cluster" not in df.columns:
        raise ValueError(f"No cluster column in {path}")
    if "sample_name" not in df.columns:
        raise ValueError(f"No sample_name in {path}")
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--clinical", type=Path, required=True)
    ap.add_argument("--cancer", "-c", type=str, required=True)
    ap.add_argument(
        "--labels",
        type=str,
        nargs="+",
        help="Pairs method_name=path e.g. SubtypeDCC=results/BRCA.dcc",
    )
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    clin = load_clinical_table(args.clinical, args.cancer)
    rows = []
    for item in args.labels:
        if "=" not in item:
            raise ValueError("Expected METHOD=path")
        name, path_s = item.split("=", 1)
        path = Path(path_s)
        lab = load_any_labels(path)
        merged = merge_labels_clinical(lab, clin)
        if merged.empty:
            continue
        nlp = survival_neg_log10_p(merged)
        n_sig, _ = clinical_significant_count(merged)
        rows.append(
            {
                "method": name,
                "neg_log10_logrank_p": nlp,
                "n_significant_clinical": n_sig,
                "n_samples": len(merged),
            }
        )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(rows, indent=2), encoding="utf-8")
    print(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()
