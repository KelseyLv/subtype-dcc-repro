"""
Figure 4 style (KIRC): (A) t-SNE of Subtype-DCC embeddings colored by cluster;
(B) two scatter plots of prominent RNA biomarkers — discrete red (native subtype)
vs blue (others), matching the paper’s Figure 4B style (not a continuous heatmap).

Usage from repo root::

    .venv\\Scripts\\python.exe repro\\plot_figure4.py -c KIRC --out repro_out/figure4_KIRC.png
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from sklearn.manifold import TSNE

REPRO_DIR = Path(__file__).resolve().parent
REPO_ROOT = REPRO_DIR.parent
if str(REPRO_DIR) not in sys.path:
    sys.path.insert(0, str(REPRO_DIR))

from data_io import resolve_fea_root  # noqa: E402


def load_rna_samples_genes(cancer: str) -> Tuple[np.ndarray, List[str], List[str]]:
    """RNA block: rows=samples, cols=genes (fea files are genes x samples)."""
    root = resolve_fea_root(REPO_ROOT)
    path = root / cancer / "rna.fea"
    df = pd.read_csv(path, sep=",", header=0, index_col=0)
    genes = list(df.index.astype(str))
    samples = list(df.columns.astype(str))
    X = df.values.T.astype(np.float64)
    return X, samples, genes


def pick_top_biomarker_pairs(
    uniq: List[int],
    markers: List[int],
    neglogp: List[float],
    marker_genes: List[str],
    k: int = 2,
) -> List[Tuple[int, int, float, str]]:
    """Pick the k subtypes whose marker genes have the strongest −log10 p (for two B panels)."""
    keyed = [
        (int(c), int(mj), float(nl), str(g))
        for c, mj, nl, g in zip(uniq, markers, neglogp, marker_genes)
    ]
    keyed.sort(key=lambda t: -t[2])
    return keyed[: min(k, len(keyed))]


def scatter_biomarker_panel(
    ax,
    X_rna: np.ndarray,
    clus: np.ndarray,
    gene_j: int,
    gene_name: str,
    home_clus: int,
    neglogp: float,
    rng: np.random.Generator,
) -> None:
    """Strip-chart style: x = subtype label + jitter, y = expression; red = native cluster."""
    y = X_rna[:, gene_j].astype(np.float64)
    labs = clus.astype(np.float64)
    jitter = rng.uniform(-0.32, 0.32, size=len(labs))
    xj = labs + jitter
    is_home = clus == home_clus
    ax.scatter(
        xj[~is_home],
        y[~is_home],
        c="#2166AC",
        s=24,
        alpha=0.75,
        edgecolors="none",
        rasterized=True,
        label="Other subtypes",
    )
    ax.scatter(
        xj[is_home],
        y[is_home],
        c="#B2182B",
        s=30,
        alpha=0.88,
        edgecolors="none",
        rasterized=True,
        label=f"Subtype {int(home_clus)} (marker host)",
    )
    ticks = sorted(np.unique(clus))
    ax.set_xticks(ticks)
    ax.set_xticklabels([str(int(t)) for t in ticks])
    ax.set_xlabel("Subtype-DCC cluster")
    ax.set_ylabel("RNA expression (fea)")
    ax.set_title(
        f"{gene_name}\n"
        f"(−log₁₀ p ≈ {neglogp:.1f}; red = higher in subtype {int(home_clus)})",
        fontsize=10,
    )
    ax.legend(loc="best", fontsize=7, framealpha=0.92)
    ax.grid(True, alpha=0.28)


def pick_one_marker_per_cluster(
    X: np.ndarray,
    labels: np.ndarray,
    min_per_group: int = 5,
) -> Tuple[List[int], List[float]]:
    """
    For each cluster id (sorted unique), run t-test vs rest on each gene;
    pick the gene with smallest p-value among those with higher mean in-cluster.
    Returns gene indices and -log10 p-values for annotation.
    """
    uniq = sorted(np.unique(labels))
    gene_idx: List[int] = []
    neglogp: List[float] = []
    n_genes = X.shape[1]

    for c in uniq:
        in_c = labels == c
        out_c = ~in_c
        if in_c.sum() < min_per_group or out_c.sum() < min_per_group:
            gene_idx.append(0)
            neglogp.append(0.0)
            continue
        best_j = -1
        best_p = 1.0
        mean_in = X[in_c].mean(axis=0)
        mean_out = X[out_c].mean(axis=0)
        for j in range(n_genes):
            a = X[in_c, j]
            b = X[out_c, j]
            if np.nanstd(a) == 0 and np.nanstd(b) == 0:
                continue
            _, p = ttest_ind(a, b, equal_var=False)
            if not np.isfinite(p):
                continue
            if mean_in[j] <= mean_out[j]:
                continue
            if p < best_p:
                best_p = float(p)
                best_j = j
        if best_j < 0:
            for j in range(n_genes):
                a = X[in_c, j]
                b = X[out_c, j]
                if np.nanstd(a) == 0 and np.nanstd(b) == 0:
                    continue
                _, p = ttest_ind(a, b, equal_var=False)
                if np.isfinite(p) and p < best_p:
                    best_p = float(p)
                    best_j = j
        gene_idx.append(best_j if best_j >= 0 else 0)
        neglogp.append(-np.log10(best_p + 1e-300))

    return gene_idx, neglogp


def plot_figure4(
    cancer: str,
    labels_path: Path,
    embedding_path: Path,
    out_path: Path,
    perplexity: float = 30.0,
    seed: int = 42,
    dpi: int = 150,
) -> None:
    labels_df = pd.read_csv(labels_path, sep="\t")
    if "dcc" in labels_df.columns:
        labels_df = labels_df.rename(columns={"dcc": "cluster"})
    sample_order = list(labels_df["sample_name"].astype(str))
    clus = labels_df["cluster"].values

    emb = pd.read_csv(embedding_path, sep="\t", index_col=0)
    emb = emb.loc[[s for s in sample_order if s in emb.index]]
    if len(emb) != len(sample_order):
        ok = set(emb.index)
        mask = [s in ok for s in sample_order]
        sample_order = [s for s, m in zip(sample_order, mask) if m]
        clus = clus[np.array(mask)]
    X_emb = emb.loc[sample_order].values.astype(np.float64)

    X_rna, rna_samples, gene_names = load_rna_samples_genes(cancer)
    samp_to_row = {s: i for i, s in enumerate(rna_samples)}
    rows = [samp_to_row[s] for s in sample_order]
    X_rna = X_rna[rows, :]

    tsne = TSNE(
        n_components=2,
        perplexity=min(perplexity, max(5, (len(X_emb) - 1) / 3)),
        random_state=seed,
        init="random",
        learning_rate="auto",
    )
    Z = tsne.fit_transform(X_emb)

    markers, neglogp = pick_one_marker_per_cluster(X_rna, clus)
    marker_genes = [gene_names[j] for j in markers]

    uniq = sorted(np.unique(clus))
    bio_pairs = pick_top_biomarker_pairs(uniq, markers, neglogp, marker_genes, k=2)
    rng = np.random.default_rng(seed)

    fig = plt.figure(figsize=(13, 6))
    gs = fig.add_gridspec(
        2,
        2,
        width_ratios=[1.28, 1.0],
        height_ratios=[1.0, 1.0],
        wspace=0.28,
        hspace=0.42,
    )
    ax_a = fig.add_subplot(gs[:, 0])
    ax_b1 = fig.add_subplot(gs[0, 1])
    ax_b2 = fig.add_subplot(gs[1, 1])

    ax = ax_a
    tab10 = plt.colormaps["tab10"].resampled(10)
    palette = [tab10(k / 9.0) for k in range(10)]
    for i, c in enumerate(uniq):
        m = clus == c
        ax.scatter(
            Z[m, 0],
            Z[m, 1],
            s=18,
            alpha=0.85,
            color=palette[i % 10],
            edgecolors="none",
            label=f"Subtype {int(c)}",
        )
    ax.set_title(f"{cancer}: t-SNE of Subtype-DCC embedding (Fig. 4A style)")
    ax.set_xlabel("t-SNE 1")
    ax.set_ylabel("t-SNE 2")
    ax.legend(loc="best", fontsize=8, framealpha=0.9)
    ax.grid(True, alpha=0.25)

    # Figure 4B: two discrete scatter plots (paper: red in native subtype, blue elsewhere)
    if len(bio_pairs) >= 1:
        c0, j0, nl0, g0 = bio_pairs[0]
        scatter_biomarker_panel(ax_b1, X_rna, clus, j0, g0, c0, nl0, rng)
        ax_b1.text(
            0.02,
            0.98,
            "B (i)",
            transform=ax_b1.transAxes,
            fontsize=11,
            fontweight="bold",
            va="top",
        )
    if len(bio_pairs) >= 2:
        c1, j1, nl1, g1 = bio_pairs[1]
        scatter_biomarker_panel(ax_b2, X_rna, clus, j1, g1, c1, nl1, rng)
        ax_b2.text(
            0.02,
            0.98,
            "B (ii)",
            transform=ax_b2.transAxes,
            fontsize=11,
            fontweight="bold",
            va="top",
        )
    else:
        ax_b2.text(0.5, 0.5, "B (ii): need ≥2 clusters", ha="center", va="center")
        ax_b2.axis("off")

    ax_a.text(
        0.02,
        0.98,
        "A",
        transform=ax_a.transAxes,
        fontsize=11,
        fontweight="bold",
        va="top",
    )

    fig.suptitle(f"{cancer}: Subtype-DCC Figure 4–style (A: embedding; B: RNA biomarkers)", y=1.02)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print("Saved", out_path)


def main() -> None:
    ap = argparse.ArgumentParser(description="Paper Figure 4–style panel for Subtype-DCC.")
    ap.add_argument("--cancer", "-c", type=str, default="KIRC")
    ap.add_argument(
        "--labels",
        type=Path,
        default=None,
        help="Default: vendor/Subtype-DCC/Subtype-DCC/results/<cancer>.dcc",
    )
    ap.add_argument(
        "--embedding",
        type=Path,
        default=None,
        help="Default: vendor/Subtype-DCC/Subtype-DCC/results/<cancer>.fea",
    )
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--perplexity", type=float, default=30.0)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    c = args.cancer.upper()
    labels_path = args.labels or (
        REPO_ROOT / "vendor" / "Subtype-DCC" / "Subtype-DCC" / "results" / f"{c}.dcc"
    )
    embedding_path = args.embedding or (
        REPO_ROOT / "vendor" / "Subtype-DCC" / "Subtype-DCC" / "results" / f"{c}.fea"
    )
    out_path = args.out or (REPO_ROOT / "repro_out" / f"figure4_{c}.png")

    plot_figure4(
        cancer=c,
        labels_path=labels_path,
        embedding_path=embedding_path,
        out_path=out_path,
        perplexity=args.perplexity,
        seed=args.seed,
    )


if __name__ == "__main__":
    main()
