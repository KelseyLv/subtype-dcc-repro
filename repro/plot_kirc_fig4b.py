"""Plot Fig 4B from cache only (no Enrichr)."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

REPRO_DIR = Path(__file__).resolve().parent
if str(REPRO_DIR) not in sys.path:
    sys.path.insert(0, str(REPRO_DIR))

from kirc_manuscript_cache_io import default_cache_dir, load_rna_marker_bundle, load_tsne_bundle  # noqa: E402
from kirc_manuscript_lib import REPO_ROOT, fig4b_tsne_marker_grid  # noqa: E402


def plot(cache_dir: Path | None = None, out_dir: Path | None = None, dpi: int = 200) -> Path:
    cache_dir = cache_dir or default_cache_dir()
    out_dir = out_dir or (REPO_ROOT / "repro_out/kirc_manuscript_figs")
    Z, labels, xlim, ylim = load_tsne_bundle(cache_dir)
    X_rna, genes, marker_idx, uniq_clusters, labels2 = load_rna_marker_bundle(cache_dir)
    if not np.array_equal(labels, labels2):
        raise ValueError("labels.npy does not match RNA bundle order")
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / "KIRC_fig4B_marker_heatmap.png"
    fig4b_tsne_marker_grid(
        Z,
        X_rna,
        labels,
        genes,
        marker_idx,
        uniq_clusters,
        out,
        xlim=xlim,
        ylim=ylim,
        dpi=dpi,
    )
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache-dir", type=Path, default=None)
    ap.add_argument("--out-dir", type=Path, default=REPO_ROOT / "repro_out/kirc_manuscript_figs")
    ap.add_argument("--dpi", type=int, default=200)
    args = ap.parse_args()
    out = plot(cache_dir=args.cache_dir, out_dir=args.out_dir, dpi=args.dpi)
    print("Wrote", out)


if __name__ == "__main__":
    main()
