"""Kaplan–Meier curves colored by cluster (Figure 3 style)."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test

from data_io import load_clinical_table
from evaluate import load_labels_dcc, merge_labels_clinical


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--labels", type=Path, required=True)
    ap.add_argument("--clinical", type=Path, required=True)
    ap.add_argument("--cancer", "-c", type=str, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    clin = load_clinical_table(args.clinical, args.cancer)
    labels = load_labels_dcc(args.labels)
    merged = merge_labels_clinical(labels, clin)
    merged = merged.dropna(subset=["days", "status", "cluster"])
    merged = merged[merged["days"] > 0]

    plt.figure(figsize=(8, 6))
    kmf = KaplanMeierFitter()
    for cl in sorted(merged["cluster"].unique()):
        m = merged["cluster"] == cl
        kmf.fit(merged.loc[m, "days"], merged.loc[m, "status"], label=f"Cluster {int(cl)}")
        kmf.plot_survival_function(ci_show=True)

    try:
        r = multivariate_logrank_test(
            merged["days"], merged["cluster"], merged["status"]
        )
        plt.title(f"{args.cancer} survival (log-rank p={r.p_value:.2e})")
    except Exception:
        plt.title(f"{args.cancer} survival")

    plt.xlabel("Time (days)")
    plt.ylabel("Survival probability")
    plt.grid(True, alpha=0.3)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(args.out, dpi=150)
    plt.close()
    print("Saved", args.out)


if __name__ == "__main__":
    main()
