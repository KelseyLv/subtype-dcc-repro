"""
One-shot heavy computation for KIRC manuscript figures: t-SNE, DE lists, Enrichr (GO/KEGG), Fig7 matrix.
Writes `repro_out/kirc_manuscript_cache/` (override with --cache-dir).

Then run individual plot scripts (e.g. `plot_kirc_fig5.py`) without re-querying Enrichr.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from sklearn.manifold import TSNE

REPRO_DIR = Path(__file__).resolve().parent
REPO_ROOT = REPRO_DIR.parent
if str(REPRO_DIR) not in sys.path:
    sys.path.insert(0, str(REPRO_DIR))

from kirc_manuscript_cache_io import default_cache_dir, save_cache_bundle  # noqa: E402
from kirc_manuscript_lib import (  # noqa: E402
    enrichr_go_combined_three_ontologies,
    load_labels,
    load_omic_transpose,
    pathway_scores,
    pick_marker_gene_per_cluster,
    top_de_gene_symbols,
    tsne_axis_limits,
)


def run_compute(
    cancer: str = "KIRC",
    labels_path: Path | None = None,
    embedding_path: Path | None = None,
    cache_dir: Path | None = None,
    seed: int = 42,
    tsne_perplexity: float = 30.0,
) -> Path:
    """Run full computation and return cache directory path."""
    cancer = cancer.upper()
    cache_dir = cache_dir or default_cache_dir()
    labels_path = labels_path or (REPO_ROOT / "vendor/Subtype-DCC/Subtype-DCC/results/KIRC.dcc")
    embedding_path = embedding_path or (REPO_ROOT / "vendor/Subtype-DCC/Subtype-DCC/results/KIRC.fea")
    cache_dir.mkdir(parents=True, exist_ok=True)

    lab_df = load_labels(labels_path)
    samples = list(lab_df["sample_name"].astype(str))
    labels = lab_df["cluster"].values.astype(int)

    emb = pd.read_csv(embedding_path, sep="\t", index_col=0)
    emb = emb.loc[[s for s in samples if s in emb.index]]
    common = list(emb.index)
    li = np.array([samples.index(s) for s in common])
    labels = labels[li]
    samples = common
    X_emb = emb.values.astype(np.float64)

    perp = min(tsne_perplexity, max(5.0, (len(X_emb) - 1) / 3.0))
    tsne = TSNE(
        n_components=2,
        perplexity=perp,
        random_state=seed,
        init="random",
        learning_rate="auto",
    )
    Z = tsne.fit_transform(X_emb)
    xlim, ylim = tsne_axis_limits(Z)

    X_rna, genes = load_omic_transpose(cancer, "rna.fea", samples)
    assert X_rna.shape[0] == len(samples)
    uniq_clusters = sorted(np.unique(labels).astype(int).tolist())
    marker_idx = pick_marker_gene_per_cluster(X_rna, labels)

    X_mir, mir_names = load_omic_transpose(cancer, "miRNA.fea", samples)

    go_parts: List[pd.DataFrame] = []
    kegg_parts: List[pd.DataFrame] = []

    try:
        import gseapy as gp
    except ImportError:
        print("gseapy not installed; GO/KEGG/Fig7 matrix will be empty.")
        gp = None  # type: ignore[assignment]

    if gp is not None:
        for c in uniq_clusters:
            glist = top_de_gene_symbols(X_rna, labels, genes, int(c), top_n=150)
            if len(glist) < 10:
                continue
            try:
                res = enrichr_go_combined_three_ontologies(glist, gp, top_terms=15)
                if not res.empty:
                    r2 = res.copy()
                    r2["cluster"] = int(c)
                    go_parts.append(r2)
            except Exception as e:  # noqa: BLE001
                print(f"GO enrich cluster {c} failed: {e}")
            try:
                enr = gp.enrichr(
                    gene_list=glist,
                    gene_sets=["KEGG_2021_Human"],
                    organism="human",
                    outdir=None,
                    verbose=False,
                )
                sub = enr.results.head(8).copy()
                sub["cluster"] = f"C{int(c)}"
                kegg_parts.append(sub)
            except Exception as e:  # noqa: BLE001
                print(f"KEGG enrich cluster {c} failed: {e}")

    go_fig5 = pd.concat(go_parts, ignore_index=True) if go_parts else pd.DataFrame()
    kegg_fig6 = pd.concat(kegg_parts, ignore_index=True) if kegg_parts else pd.DataFrame()

    pathway_terms: List[str] = []
    if not kegg_fig6.empty:
        pathway_terms = list(dict.fromkeys(kegg_fig6["Term"].tolist()))[:12]

    var = X_mir.var(axis=0)
    top_m = np.argsort(-var)[: min(40, X_mir.shape[1])]
    mir_sel = [mir_names[i] for i in top_m]

    lib_kegg = None
    if gp is not None:
        for org in ("Human", "human"):
            try:
                lib_kegg = gp.get_library(name="KEGG_2021_Human", organism=org)
                break
            except Exception:
                continue

    mat = np.zeros((len(mir_sel), len(pathway_terms)))
    mat[:] = np.nan
    if gp is not None and lib_kegg is not None and pathway_terms:
        for pi, term in enumerate(pathway_terms):
            pgenes = lib_kegg.get(term, [])
            ps = pathway_scores(X_rna, samples, genes, pgenes)
            if np.all(np.isnan(ps)):
                continue
            for mi, mir_i in enumerate(top_m):
                v = X_mir[:, mir_i]
                m = np.isfinite(ps) & np.isfinite(v)
                if m.sum() < 15:
                    continue
                r, p = pearsonr(v[m], ps[m])
                mat[mi, pi] = -np.log10(p + 1e-300) * np.sign(r) if np.isfinite(p) else 0.0

    meta = {
        "cancer": cancer,
        "labels_path": str(labels_path),
        "embedding_path": str(embedding_path),
        "tsne_seed": seed,
        "tsne_perplexity": float(perp),
        "n_samples": len(samples),
    }

    save_cache_bundle(
        cache_dir=cache_dir,
        meta=meta,
        Z=Z,
        labels=labels,
        samples=samples,
        xlim=xlim,
        ylim=ylim,
        X_rna=X_rna,
        genes=genes,
        marker_idx=marker_idx,
        uniq_clusters=uniq_clusters,
        go_fig5=go_fig5,
        kegg_fig6=kegg_fig6,
        fig7_matrix=mat,
        fig7_pathways=pathway_terms,
        fig7_mirlabels=mir_sel,
        X_mir=X_mir,
        mir_names=mir_names,
    )
    (cache_dir / "README.txt").write_text(
        "KIRC manuscript cache. Regenerate with:\n"
        f"  python repro/kirc_manuscript_compute_cache.py --cache-dir {cache_dir}\n"
        "Plot Fig 5 only from cache:\n"
        f"  python repro/plot_kirc_fig5.py --cache-dir {cache_dir}\n"
        "Cache tables: go_fig5.csv, kegg_fig6.csv (no pyarrow required).\n",
        encoding="utf-8",
    )
    print("Saved cache to", cache_dir)
    return cache_dir


def main() -> None:
    ap = argparse.ArgumentParser(description="Compute and save KIRC manuscript figure cache.")
    ap.add_argument("--cancer", type=str, default="KIRC")
    ap.add_argument(
        "--labels",
        type=Path,
        default=REPO_ROOT / "vendor/Subtype-DCC/Subtype-DCC/results/KIRC.dcc",
    )
    ap.add_argument(
        "--embedding",
        type=Path,
        default=REPO_ROOT / "vendor/Subtype-DCC/Subtype-DCC/results/KIRC.fea",
    )
    ap.add_argument("--cache-dir", type=Path, default=None, help="Default: repro_out/kirc_manuscript_cache")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--tsne-perplexity", type=float, default=30.0)
    args = ap.parse_args()
    run_compute(
        cancer=args.cancer,
        labels_path=args.labels,
        embedding_path=args.embedding,
        cache_dir=args.cache_dir,
        seed=args.seed,
        tsne_perplexity=args.tsne_perplexity,
    )


if __name__ == "__main__":
    main()
