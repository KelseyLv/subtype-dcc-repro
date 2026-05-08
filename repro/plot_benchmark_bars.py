"""Bar chart from run_benchmark_eval.py JSON (performance comparison)."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary-json", type=Path, required=True)
    ap.add_argument("--out-prefix", type=Path, required=True)
    args = ap.parse_args()

    rows = json.loads(args.summary_json.read_text(encoding="utf-8"))
    methods = [r["method"] for r in rows]
    x = np.arange(len(methods))

    fig, ax = plt.subplots(figsize=(max(8, len(methods) * 0.5), 5))
    ax.bar(x, [r["neg_log10_logrank_p"] for r in rows], color="steelblue")
    ax.set_xticks(x)
    ax.set_xticklabels(methods, rotation=45, ha="right")
    ax.set_ylabel("-log10 log-rank p")
    ax.set_title("Survival separation")
    plt.tight_layout()
    p1 = Path(str(args.out_prefix) + "_survival.png")
    p1.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(p1, dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(max(8, len(methods) * 0.5), 5))
    ax.bar(x, [r["n_significant_clinical"] for r in rows], color="darkorange")
    ax.set_xticks(x)
    ax.set_xticklabels(methods, rotation=45, ha="right")
    ax.set_ylabel("# significant clinical variables (of 6)")
    ax.set_title("Clinical enrichment")
    plt.tight_layout()
    p2 = Path(str(args.out_prefix) + "_clinical.png")
    fig.savefig(p2, dpi=150)
    plt.close(fig)
    print("Saved", p1, p2)


if __name__ == "__main__":
    main()
