"""
Figure 3 style: 3×3 Kaplan–Meier survival curves by Subtype-DCC cluster for nine TCGA cancers.

Uses the same clinical merge and multivariate log-rank test as ``plot_survival.py`` /
``evaluate.py``. Subplot titles are cancer acronyms; log-rank *p* is annotated on each panel.

Usage from repository root::

    .venv\\Scripts\\python.exe repro\\plot_figure3_survival_grid.py --out repro_out/figure3_survival_9cancers.png
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List, Sequence

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test

REPRO_DIR = Path(__file__).resolve().parent
REPO_ROOT = REPRO_DIR.parent
if str(REPRO_DIR) not in sys.path:
    sys.path.insert(0, str(REPRO_DIR))

from data_io import load_clinical_table  # noqa: E402
from evaluate import load_labels_dcc, merge_labels_clinical  # noqa: E402

# Colorblind-friendly palette (Wong / similar); up to 5 clusters in Table 1 setting.
CLUSTER_LINE_COLORS: Sequence[str] = (
    "#0072B2",  # blue
    "#D55E00",  # vermillion
    "#009E73",  # bluish green
    "#CC79A7",  # reddish purple
    "#E69F00",  # orange
)

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


def _prepare_merged(
    cancer: str,
    labels_path: Path,
    clinical_path: Path,
) -> pd.DataFrame:
    clin = load_clinical_table(clinical_path, cancer)
    labels = load_labels_dcc(labels_path)
    merged = merge_labels_clinical(labels, clin)
    merged = merged.dropna(subset=["days", "status", "cluster"])
    merged = merged[merged["days"] > 0]
    return merged


def _format_p(p: float) -> str:
    if not np.isfinite(p) or p <= 0:
        return r"$p < 10^{-300}$"
    if p < 1e-4:
        return rf"$p = {p:.2e}$"
    return rf"$p = {p:.4g}$"


def plot_one_cancer_ax(
    ax: mpl.axes.Axes,
    cancer: str,
    labels_path: Path,
    clinical_path: Path,
    line_colors: Sequence[str],
    show_legend: bool,
    show_xlabel: bool,
    show_ylabel: bool,
) -> None:
    merged = _prepare_merged(cancer, labels_path, clinical_path)
    ax.set_title(cancer, fontsize=12, fontweight="semibold")
    if merged.empty or merged["cluster"].nunique() < 2:
        ax.text(0.5, 0.5, "Insufficient overlap\nor single cluster", ha="center", va="center", transform=ax.transAxes)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")
        return

    clusters = sorted(merged["cluster"].unique())
    for j, cl in enumerate(clusters):
        m = merged["cluster"] == cl
        color = line_colors[j % len(line_colors)]
        kmf = KaplanMeierFitter()
        kmf.fit(
            merged.loc[m, "days"],
            merged.loc[m, "status"],
            label=f"Cluster {int(cl)}",
        )
        kmf.plot_survival_function(
            ax=ax,
            ci_show=True,
            color=color,
            show_censors=False,
            linewidth=2.0,
        )

    try:
        r = multivariate_logrank_test(
            merged["days"],
            merged["cluster"],
            merged["status"],
        )
        p_val = float(r.p_value)
        p_str = _format_p(p_val)
    except Exception:
        p_str = r"$p$: n/a"

    ax.text(
        0.98,
        0.05,
        f"log-rank\n{p_str}",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=9,
        linespacing=1.15,
        bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="#cccccc", alpha=0.92),
    )
    ax.grid(True, alpha=0.35)
    ax.set_xlim(left=0)
    ax.set_ylim(0, 1.02)
    if show_xlabel:
        ax.set_xlabel("Time (days)", fontsize=10)
    else:
        ax.set_xlabel("")
    if show_ylabel:
        ax.set_ylabel("Survival probability", fontsize=10)
    else:
        ax.set_ylabel("")
    if show_legend:
        ax.legend(loc="upper right", fontsize=7, framealpha=0.9)


def main() -> None:
    ap = argparse.ArgumentParser(description="3×3 KM curves (Figure 3 style, nine cancers).")
    ap.add_argument(
        "--clinical",
        type=Path,
        default=REPO_ROOT / "data" / "clinical" / "cBioportal_data.tsv",
    )
    ap.add_argument(
        "--labels-dir",
        type=Path,
        default=REPO_ROOT / "vendor" / "Subtype-DCC" / "Subtype-DCC" / "results",
        help="Directory containing <CANCER>.dcc TSV files.",
    )
    ap.add_argument(
        "--cancers",
        nargs="*",
        default=DEFAULT_CANCERS,
        help="Exactly nine cancer codes in row-major order for the 3×3 grid (default: Table 1 set).",
    )
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--dpi", type=int, default=200)
    ap.add_argument("--fig-width", type=float, default=14.0)
    ap.add_argument("--fig-height", type=float, default=12.0)
    args = ap.parse_args()

    cancers: List[str] = list(args.cancers)
    if len(cancers) != 9:
        raise SystemExit("Please pass exactly nine cancers for the 3×3 layout.")

    # Unified typography
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.titlesize": 12,
            "axes.labelsize": 10,
            "legend.fontsize": 8,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )

    fig, axes = plt.subplots(3, 3, figsize=(args.fig_width, args.fig_height), sharex=False, sharey=True)
    axes_flat = axes.flatten()

    for i, cancer in enumerate(cancers):
        ax = axes_flat[i]
        labels_path = args.labels_dir / f"{cancer}.dcc"
        if not labels_path.is_file():
            ax.set_title(cancer, fontsize=12, fontweight="semibold")
            ax.text(0.5, 0.5, f"Missing\n{labels_path.name}", ha="center", va="center", transform=ax.transAxes)
            ax.axis("off")
            continue
        row, col = divmod(i, 3)
        show_xlabel = row == 2
        show_ylabel = col == 0
        show_legend = i == 0
        plot_one_cancer_ax(
            ax,
            cancer,
            labels_path,
            args.clinical,
            CLUSTER_LINE_COLORS,
            show_legend=show_legend,
            show_xlabel=show_xlabel,
            show_ylabel=show_ylabel,
        )

    fig.suptitle(
        "Kaplan–Meier overall survival by Subtype-DCC cluster (TCGA)",
        fontsize=14,
        fontweight="semibold",
        y=1.02,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=args.dpi, bbox_inches="tight")
    plt.close(fig)
    print("Saved", args.out.resolve())


if __name__ == "__main__":
    main()
