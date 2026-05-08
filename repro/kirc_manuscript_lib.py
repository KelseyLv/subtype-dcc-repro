"""
Shared loaders, enrichment helpers, and matplotlib figure builders for KIRC manuscript-style plots (Fig 4–7).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
import seaborn as sns
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.stats import pearsonr, ttest_ind
from sklearn.manifold import TSNE

REPRO_DIR = Path(__file__).resolve().parent
REPO_ROOT = REPRO_DIR.parent
if str(REPRO_DIR) not in sys.path:
    sys.path.insert(0, str(REPRO_DIR))

from data_io import resolve_fea_root  # noqa: E402

CACHE_VERSION = 1


def load_labels(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    if "dcc" in df.columns:
        df = df.rename(columns={"dcc": "cluster"})
    return df


def load_omic_transpose(cancer: str, fname: str, samples_ref: List[str]) -> Tuple[np.ndarray, List[str]]:
    root = resolve_fea_root(REPO_ROOT)
    p = root / cancer / fname
    df = pd.read_csv(p, sep=",", header=0, index_col=0)
    common = [s for s in samples_ref if s in df.columns]
    df = df[common]
    feats = [str(x) for x in df.index]
    return df.values.T.astype(np.float64), feats


def pick_marker_gene_per_cluster(X: np.ndarray, labels: np.ndarray, min_n: int = 8) -> List[int]:
    uniq = sorted(np.unique(labels))
    picks: List[int] = []
    for c in uniq:
        in_c = labels == c
        out_c = ~in_c
        if in_c.sum() < min_n or out_c.sum() < min_n:
            picks.append(0)
            continue
        best_j, best_p = -1, 1.0
        mean_in = X[in_c].mean(axis=0)
        mean_out = X[out_c].mean(axis=0)
        for j in range(X.shape[1]):
            if mean_in[j] <= mean_out[j]:
                continue
            a, b = X[in_c, j], X[out_c, j]
            if np.nanstd(a) < 1e-12 and np.nanstd(b) < 1e-12:
                continue
            _, p = ttest_ind(a, b, equal_var=False)
            if np.isfinite(p) and p < best_p:
                best_p, best_j = float(p), j
        if best_j < 0:
            for j in range(X.shape[1]):
                a, b = X[in_c, j], X[out_c, j]
                if np.nanstd(a) < 1e-12 and np.nanstd(b) < 1e-12:
                    continue
                _, p = ttest_ind(a, b, equal_var=False)
                if np.isfinite(p) and p < best_p:
                    best_p, best_j = float(p), j
        picks.append(best_j if best_j >= 0 else 0)
    return picks


def tsne_axis_limits(Z: np.ndarray, pad_frac: float = 0.05) -> Tuple[Tuple[float, float], Tuple[float, float]]:
    x0, x1 = float(Z[:, 0].min()), float(Z[:, 0].max())
    y0, y1 = float(Z[:, 1].min()), float(Z[:, 1].max())
    dx = max(x1 - x0, 1e-6)
    dy = max(y1 - y0, 1e-6)
    px = pad_frac * dx
    py = pad_frac * dy
    return (x0 - px, x1 + px), (y0 - py, y1 + py)


def overlap_count_from_overlap(s: str) -> int:
    s = str(s)
    if "/" not in s:
        return 0
    a, _ = s.split("/", 1)
    try:
        return int(a.strip())
    except ValueError:
        return 0


GO_ENRICHR_LIBS: Tuple[Tuple[str, str], ...] = (
    ("GO_Biological_Process_2021", "BP"),
    ("GO_Cellular_Component_2021", "CC"),
    ("GO_Molecular_Function_2021", "MF"),
)

ONTOLOGY_COLORS = {"BP": "#2171b5", "CC": "#31a354", "MF": "#fd8d3c"}


def enrichr_go_combined_three_ontologies(
    glist: List[str], gp, top_terms: int = 15
) -> pd.DataFrame:
    parts: List[pd.DataFrame] = []
    for lib, onto in GO_ENRICHR_LIBS:
        try:
            enr = gp.enrichr(
                gene_list=glist,
                gene_sets=[lib],
                organism="human",
                outdir=None,
                verbose=False,
            )
            r = enr.results.copy()
            r["ontology"] = onto
            parts.append(r)
        except Exception:
            continue
    if not parts:
        return pd.DataFrame()
    comb = pd.concat(parts, ignore_index=True)
    if comb.empty or "Adjusted P-value" not in comb.columns:
        return pd.DataFrame()
    comb["gene_count"] = comb["Overlap"].astype(str).map(overlap_count_from_overlap)
    comb["neglogpadj"] = -np.log10(comb["Adjusted P-value"].astype(float) + 1e-300)
    comb = comb.sort_values("neglogpadj", ascending=False).drop_duplicates(subset=["Term"], keep="first")
    comb = comb.nlargest(top_terms, "neglogpadj")
    comb = comb.sort_values("neglogpadj", ascending=True)
    comb["color"] = comb["ontology"].map(lambda o: ONTOLOGY_COLORS.get(str(o), "#999999"))
    return comb


def shorten_go_term(t: str, max_len: int = 52) -> str:
    t = str(t).strip()
    if len(t) <= max_len:
        return t
    return t[: max_len - 1] + "…"


def top_de_gene_symbols(
    X: np.ndarray, labels: np.ndarray, genes: List[str], cluster: int, top_n: int = 120
) -> List[str]:
    in_c = labels == cluster
    out_c = ~in_c
    scores = []
    for j in range(X.shape[1]):
        a, b = X[in_c, j], X[out_c, j]
        if np.nanstd(a) < 1e-12:
            continue
        t, p = ttest_ind(a, b, equal_var=False)
        if not np.isfinite(p):
            continue
        if X[in_c, j].mean() <= X[out_c, j].mean():
            continue
        scores.append((float(t), -np.log10(p + 1e-300), j))
    scores.sort(key=lambda x: (-x[0], -x[1]))
    out: List[str] = []
    for _, __, j in scores[:top_n]:
        sym = str(genes[j]).split("|")[0].strip()
        if sym and sym not in out:
            out.append(sym)
    return out[:top_n]


def pathway_scores(
    X_rna: np.ndarray, samples: List[str], gene_names: List[str], pathway_genes: List[str]
) -> np.ndarray:
    sym_to_j = {str(g).split("|")[0].upper(): j for j, g in enumerate(gene_names)}
    js = [sym_to_j[g.upper()] for g in pathway_genes if g.upper() in sym_to_j]
    if len(js) < 3:
        return np.full(X_rna.shape[0], np.nan)
    sub = X_rna[:, js]
    z = (sub - sub.mean(axis=0, keepdims=True)) / (sub.std(axis=0, keepdims=True) + 1e-8)
    return np.nanmean(z, axis=1)


def fig4a_tsne(
    Z: np.ndarray,
    labels: np.ndarray,
    out: Path,
    dpi: int = 200,
    xlim: Optional[Tuple[float, float]] = None,
    ylim: Optional[Tuple[float, float]] = None,
) -> None:
    uniq = sorted(np.unique(labels))
    fig, ax = plt.subplots(figsize=(7.2, 6.2), facecolor="white")
    ax.set_facecolor("#fafafa")
    cmap = plt.colormaps["Set1"].resampled(max(len(uniq), 3))
    for i, c in enumerate(uniq):
        m = labels == c
        ax.scatter(
            Z[m, 0],
            Z[m, 1],
            s=36,
            alpha=0.88,
            c=[cmap(i % 9)] * int(m.sum()),
            edgecolors="white",
            linewidths=0.35,
            label=f"Cluster {int(c)}",
        )
    ax.set_xlabel("t-SNE 1", fontsize=12)
    ax.set_ylabel("t-SNE 2", fontsize=12)
    ax.set_title("KIRC — Subtype-DCC embedding (t-SNE)", fontsize=13, fontweight="bold")
    ax.legend(frameon=True, fancybox=False, shadow=False, loc="best")
    ax.grid(False)
    sns.despine(ax=ax)
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def fig4b_tsne_marker_grid(
    Z: np.ndarray,
    X_rna: np.ndarray,
    labels: np.ndarray,
    gene_names: List[str],
    marker_idx: List[int],
    uniq_clusters: List[int],
    out: Path,
    xlim: Tuple[float, float],
    ylim: Tuple[float, float],
    expr_vmin: float = 0.0,
    expr_vmax: float = 1.5,
    dpi: int = 200,
) -> None:
    mpl_markers = ["o", "^", "s", "+"]
    marker_by_cluster = {int(c): mpl_markers[i % 4] for i, c in enumerate(uniq_clusters)}
    cmap = plt.get_cmap("coolwarm")
    norm = mcolors.Normalize(vmin=expr_vmin, vmax=expr_vmax)
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    fig, axes = plt.subplots(2, 2, figsize=(12.5, 11), facecolor="white", layout="constrained")
    axes_flat = axes.ravel()

    for ax_i, ax in enumerate(axes_flat):
        if ax_i >= len(uniq_clusters):
            ax.axis("off")
            continue
        c_targ = int(uniq_clusters[ax_i])
        gj = marker_idx[ax_i]
        expr = np.clip(X_rna[:, gj], expr_vmin, expr_vmax)
        sym = str(gene_names[gj]).split("|")[0].strip() or f"gene{gj}"
        draw_order = [c for c in uniq_clusters if int(c) != c_targ] + [c_targ]
        for c_cl in draw_order:
            c_int = int(c_cl)
            m = labels == c_int
            if not np.any(m):
                continue
            mk = marker_by_cluster.get(c_int, "o")
            s_ = 44 if c_int == c_targ else 36
            if mk == "+":
                ax.scatter(
                    Z[m, 0],
                    Z[m, 1],
                    c=expr[m],
                    cmap=cmap,
                    norm=norm,
                    marker=mk,
                    s=s_,
                    linewidths=0.9,
                    alpha=0.9,
                )
            else:
                ax.scatter(
                    Z[m, 0],
                    Z[m, 1],
                    c=expr[m],
                    cmap=cmap,
                    norm=norm,
                    marker=mk,
                    s=s_,
                    edgecolors="white",
                    linewidths=0.3,
                    alpha=0.9,
                )
        ax.set_xlabel("t-SNE 1", fontsize=11)
        ax.set_ylabel("t-SNE 2", fontsize=11)
        ax.set_title(f"{sym}(cluster{c_targ})", fontsize=12, fontweight="bold")
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_facecolor("#fafafa")
        ax.grid(False)
        sns.despine(ax=ax)
        cax = inset_axes(ax, width="36%", height="3.8%", loc="upper left", borderpad=0.8)
        cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
        cb.ax.tick_params(labelsize=7, length=2, pad=1)
        cb.set_ticks([expr_vmin, 0.5, 1.0, expr_vmax])

    fig.suptitle("KIRC — Fig. 4B: subtype markers on t-SNE (same embedding as 4A)", fontsize=13, fontweight="bold")
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def fig5_go_bars_from_table(go_fig5: pd.DataFrame, uniq_clusters: List[int], out: Path, dpi: int = 200) -> None:
    """go_fig5: columns include cluster, Term, gene_count, color (optional), ontology."""
    fig5, axes5 = plt.subplots(2, 2, figsize=(14, 12), facecolor="white", layout="constrained")
    axes5_flat = axes5.ravel()
    legend_handles = [
        Patch(facecolor=ONTOLOGY_COLORS["BP"], edgecolor="none", label="BP"),
        Patch(facecolor=ONTOLOGY_COLORS["CC"], edgecolor="none", label="CC"),
        Patch(facecolor=ONTOLOGY_COLORS["MF"], edgecolor="none", label="MF"),
    ]
    for ax_i, c in enumerate(uniq_clusters):
        ax = axes5_flat[ax_i]
        ax.set_facecolor("white")
        if go_fig5.empty or "cluster" not in go_fig5.columns:
            ax.text(0.5, 0.5, "No GO cache", ha="center", va="center", transform=ax.transAxes)
            ax.axis("off")
            continue
        res = go_fig5[go_fig5["cluster"] == int(c)].copy()
        if res.empty:
            ax.text(0.5, 0.5, f"Cluster {c}: no rows", ha="center", va="center", transform=ax.transAxes)
            ax.axis("off")
            continue
        if "color" not in res.columns:
            res["color"] = res["ontology"].map(lambda o: ONTOLOGY_COLORS.get(str(o), "#999999"))
        y_pos = np.arange(len(res))
        gc = res["gene_count"].astype(float).values
        ax.barh(y_pos, gc, height=0.62, color=res["color"].values, edgecolor="none", linewidth=0)
        ax.set_yticks(y_pos)
        ax.set_yticklabels([shorten_go_term(t) for t in res["Term"].values], fontsize=7)
        ax.set_xlabel("Gene Count", fontsize=10)
        ax.set_ylabel("")
        ax.set_title(f"GO enrichment — cluster {c}", fontsize=11, fontweight="bold", pad=6)
        xmax = max(float(gc.max()), 1.0) * 1.32
        ax.set_xlim(0, xmax)
        ax.tick_params(axis="x", labelsize=9)
        ax.grid(False)
        sns.despine(ax=ax)
        ax.legend(
            handles=legend_handles,
            loc="upper right",
            frameon=True,
            fancybox=False,
            edgecolor="#cccccc",
            fontsize=8,
            title="GO",
            title_fontsize=8,
            borderpad=0.35,
        )
    fig5.suptitle(
        "KIRC — Fig. 5: GO enrichment by subtype (Biological Process / Cellular Component / Molecular Function)",
        fontsize=13,
        fontweight="bold",
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    fig5.savefig(out, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig5)


def fig6_kegg_dot(kegg_fig6: pd.DataFrame, out: Path, dpi: int = 200) -> None:
    fig6, ax6 = plt.subplots(figsize=(10, 7), facecolor="white")
    if not kegg_fig6.empty:
        comb = kegg_fig6.copy()
        comb["neglogp"] = -np.log10(comb["Adjusted P-value"].astype(float) + 1e-300)
        comb["GeneRatio"] = comb["Overlap"].astype(str).apply(
            lambda s: float(s.split("/")[0]) / max(float(s.split("/")[1]), 1.0) if "/" in s else 0.05
        )
        sns.scatterplot(
            data=comb,
            x="cluster",
            y="Term",
            size="GeneRatio",
            hue="neglogp",
            palette="viridis",
            sizes=(40, 220),
            ax=ax6,
            legend="brief",
        )
        ax6.set_title("KIRC — Fig. 6 style: KEGG (top terms per cluster)", fontsize=13, fontweight="bold")
        ax6.set_xlabel("Subtype-DCC cluster")
        ax6.set_ylabel("KEGG pathway")
    else:
        ax6.text(0.5, 0.5, "No KEGG enrichment rows", ha="center")
    fig6.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig6.savefig(out, dpi=dpi, bbox_inches="tight")
    plt.close(fig6)


def fig7_mirna_pathway_heatmap(
    mat: np.ndarray,
    pathway_terms: List[str],
    mir_rowlabels: List[str],
    out: Path,
    dpi: int = 200,
) -> None:
    fig7, ax7 = plt.subplots(figsize=(max(8, len(pathway_terms) * 0.55), max(6, len(mir_rowlabels) * 0.22)))
    if pathway_terms and np.any(np.isfinite(mat)):
        sns.heatmap(
            mat,
            xticklabels=[t[:40] + "…" if len(t) > 40 else t for t in pathway_terms],
            yticklabels=[m[:28] for m in mir_rowlabels],
            cmap="RdBu_r",
            center=0,
            ax=ax7,
            cbar_kws={"label": r"sign($r$) $\times$ $-\log_{10} p$"},
        )
        ax7.set_title(
            "KIRC — Fig. 7 style: miRNA vs KEGG pathway association\n"
            "(corr with pathway mean z-score; not DIANA-miRPath)",
            fontsize=11,
        )
        plt.setp(ax7.get_xticklabels(), rotation=45, ha="right")
    else:
        ax7.text(0.5, 0.5, "Could not build miRNA–pathway matrix", ha="center")
    fig7.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig7.savefig(out, dpi=dpi, bbox_inches="tight")
    plt.close(fig7)
