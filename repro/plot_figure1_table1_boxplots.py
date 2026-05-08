"""
Figure 1–style boxplots (paper Table 1 reproduction): five methods × nine cancers.

Reads repro_out/table1_reproduction_metrics.json (from batch_evaluate_table1_methods.py).
Writes two PNGs under repro_out/:
  - figure1_neg_log10_logrank_p_boxplot.png
  - figure1_n_significant_clinical_boxplot.png

Optional: --combined repro_out/figure1_table1_boxplots_combined.png (1×2 panels).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

REPRO_DIR = Path(__file__).resolve().parent
REPO_ROOT = REPRO_DIR.parent
if str(REPRO_DIR) not in sys.path:
    sys.path.insert(0, str(REPRO_DIR))

from batch_evaluate import DEFAULT_CANCERS  # noqa: E402

METHOD_ORDER = ["subtype_dcc", "snf", "kmeans", "spectral", "nmf"]
METHOD_LABELS = {
    "subtype_dcc": "Subtype-DCC",
    "snf": "SNF",
    "kmeans": "K-means",
    "spectral": "Spectral",
    "nmf": "NMF",
}


def load_frame(metrics_path: Path) -> pd.DataFrame:
    data = json.loads(metrics_path.read_text(encoding="utf-8"))
    rows = [r for r in data.get("rows", []) if r.get("status") == "ok"]
    if len(rows) < 9 * 5:
        raise SystemExit(
            f"Expected 45 ok rows in {metrics_path}, got {len(rows)}. Run repro/batch_evaluate_table1_methods.py first."
        )
    df = pd.DataFrame(rows)
    df = df[df["method"].isin(METHOD_ORDER)]
    df["method_label"] = df["method"].map(METHOD_LABELS)
    cat = pd.Categorical(
        df["method_label"],
        categories=[METHOD_LABELS[m] for m in METHOD_ORDER],
        ordered=True,
    )
    df["method_label"] = cat
    return df


def _style_axes(ax: plt.Axes, ylabel: str, title: str) -> None:
    ax.set_xlabel("Method")
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontweight="bold")
    ax.tick_params(axis="x", rotation=15)
    sns.despine(ax=ax, trim=True)


def plot_logrank_box(df: pd.DataFrame, out: Path, dpi: int) -> None:
    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    sns.boxplot(
        data=df,
        x="method_label",
        y="neg_log10_logrank_p",
        order=[METHOD_LABELS[m] for m in METHOD_ORDER],
        width=0.55,
        color="#a8c5e0",
        fliersize=4,
        ax=ax,
    )
    sns.stripplot(
        data=df,
        x="method_label",
        y="neg_log10_logrank_p",
        order=[METHOD_LABELS[m] for m in METHOD_ORDER],
        color="0.2",
        alpha=0.65,
        size=5,
        jitter=0.18,
        ax=ax,
    )
    _style_axes(
        ax,
        ylabel=r"$-\log_{10}(\mathrm{log\text{-}rank}\;p)$",
        title="Survival separation (nine cohorts)",
    )
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def plot_clinical_box(df: pd.DataFrame, out: Path, dpi: int) -> None:
    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    sns.boxplot(
        data=df,
        x="method_label",
        y="n_significant_clinical",
        order=[METHOD_LABELS[m] for m in METHOD_ORDER],
        width=0.55,
        color="#c5e0a8",
        fliersize=4,
        ax=ax,
    )
    sns.stripplot(
        data=df,
        x="method_label",
        y="n_significant_clinical",
        order=[METHOD_LABELS[m] for m in METHOD_ORDER],
        color="0.2",
        alpha=0.65,
        size=5,
        jitter=0.18,
        ax=ax,
    )
    ax.set_ylim(-0.35, 6.35)
    ax.set_yticks(range(0, 7))
    _style_axes(
        ax,
        ylabel="Number of significant clinical variables (of six)",
        title="Clinical association (nine cohorts)",
    )
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def plot_combined(df: pd.DataFrame, out: Path, dpi: int) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.2))
    for ax, col, title, color in zip(
        axes,
        ["neg_log10_logrank_p", "n_significant_clinical"],
        [
            r"(A) $-\log_{10}(\mathrm{log\text{-}rank}\;p)$",
            "(B) Significant clinical variables (of six)",
        ],
        ["#a8c5e0", "#c5e0a8"],
    ):
        sns.boxplot(
            data=df,
            x="method_label",
            y=col,
            order=[METHOD_LABELS[m] for m in METHOD_ORDER],
            width=0.52,
            color=color,
            fliersize=3,
            ax=ax,
        )
        sns.stripplot(
            data=df,
            x="method_label",
            y=col,
            order=[METHOD_LABELS[m] for m in METHOD_ORDER],
            color="0.2",
            alpha=0.6,
            size=4,
            jitter=0.16,
            ax=ax,
        )
        ax.set_xlabel("Method")
        ax.tick_params(axis="x", rotation=15)
        ax.set_title(title, fontweight="bold")
        sns.despine(ax=ax, trim=True)
        if col == "n_significant_clinical":
            ax.set_ylim(-0.35, 6.35)
            ax.set_yticks(range(0, 7))
            ax.set_ylabel("Count")
        else:
            ax.set_ylabel(r"$-\log_{10}(p)$")
    cohorts = ", ".join(DEFAULT_CANCERS)
    fig.suptitle(
        f"Table 1 reproduction — nine cohorts ({cohorts})",
        fontsize=10,
        y=1.02,
    )
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description="Figure 1–style boxplots for five methods × nine cancers.")
    ap.add_argument(
        "--metrics-json",
        type=Path,
        default=REPO_ROOT / "repro_out" / "table1_reproduction_metrics.json",
    )
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=REPO_ROOT / "repro_out",
        help="Directory for PNG outputs.",
    )
    ap.add_argument("--dpi", type=int, default=200)
    ap.add_argument(
        "--combined",
        action="store_true",
        help="Also write figure1_table1_boxplots_combined.png (1×2).",
    )
    args = ap.parse_args()

    df = load_frame(args.metrics_json)
    out_dir = args.out_dir
    p1 = out_dir / "figure1_neg_log10_logrank_p_boxplot.png"
    p2 = out_dir / "figure1_n_significant_clinical_boxplot.png"
    plot_logrank_box(df, p1, args.dpi)
    plot_clinical_box(df, p2, args.dpi)
    print("Wrote", p1)
    print("Wrote", p2)
    if args.combined:
        p3 = out_dir / "figure1_table1_boxplots_combined.png"
        plot_combined(df, p3, args.dpi)
        print("Wrote", p3)


if __name__ == "__main__":
    main()
