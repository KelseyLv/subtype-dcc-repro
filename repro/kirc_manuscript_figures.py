"""
KIRC manuscript-style figures (4–7).

**Recommended workflow (fast iteration on plots):**

1. Run heavy compute once (t-SNE, Enrichr, Fig7 matrix)::

    .venv\\Scripts\\python.exe repro\\kirc_manuscript_compute_cache.py --cache-dir repro_out\\kirc_manuscript_cache

2. Regenerate individual figures from cache without re-calling Enrichr::

    .venv\\Scripts\\python.exe repro\\plot_kirc_fig5.py --cache-dir repro_out\\kirc_manuscript_cache

**This entrypoint** can either mirror the old “all in one” run, or plot a subset from cache::

    # Full recompute + all figures (slow)
    .venv\\Scripts\\python.exe repro\\kirc_manuscript_figures.py

    # From cache only — e.g. Fig 5 only
    .venv\\Scripts\\python.exe repro\\kirc_manuscript_figures.py --from-cache --figures fig5
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Callable, Dict, List

REPRO_DIR = Path(__file__).resolve().parent
REPO_ROOT = REPRO_DIR.parent
if str(REPRO_DIR) not in sys.path:
    sys.path.insert(0, str(REPRO_DIR))

from kirc_manuscript_cache_io import default_cache_dir  # noqa: E402
from kirc_manuscript_compute_cache import run_compute  # noqa: E402


def _write_readme(out_dir: Path, outputs: List[str], meta: dict) -> None:
    (out_dir / "README_figures.txt").write_text(
        "KIRC manuscript-style figures.\n"
        "Fig 7 uses correlation-based proxy; paper used DIANA-miRPath.\n"
        "Prefer: repro/kirc_manuscript_compute_cache.py then repro/plot_kirc_fig*.py\n"
        f"Meta: {meta}\n"
        f"Outputs: {outputs}\n",
        encoding="utf-8",
    )


def main() -> None:
    ap = argparse.ArgumentParser(description="KIRC manuscript figures (full run or cache-backed subset).")
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
    ap.add_argument("--out-dir", type=Path, default=REPO_ROOT / "repro_out/kirc_manuscript_figs")
    ap.add_argument("--cache-dir", type=Path, default=None, help="Default: repro_out/kirc_manuscript_cache")
    ap.add_argument("--from-cache", action="store_true", help="Skip compute; read cache and plot only.")
    ap.add_argument(
        "--figures",
        type=str,
        default="fig4a,fig4b,fig5,fig6,fig7",
        help="Comma-separated: fig4a,fig4b,fig5,fig6,fig7",
    )
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--tsne-perplexity", type=float, default=30.0)
    ap.add_argument("--dpi", type=int, default=200)
    args = ap.parse_args()

    cache_dir = args.cache_dir or default_cache_dir()
    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    import plot_kirc_fig4a  # noqa: PLC0415
    import plot_kirc_fig4b
    import plot_kirc_fig5
    import plot_kirc_fig6
    import plot_kirc_fig7

    plotters: Dict[str, Callable[..., Path]] = {
        "fig4a": plot_kirc_fig4a.plot,
        "fig4b": plot_kirc_fig4b.plot,
        "fig5": plot_kirc_fig5.plot,
        "fig6": plot_kirc_fig6.plot,
        "fig7": plot_kirc_fig7.plot,
    }

    names = [x.strip().lower() for x in args.figures.split(",") if x.strip()]
    for n in names:
        if n not in plotters:
            raise SystemExit(f"Unknown figure key: {n}. Choose from {sorted(plotters)}")

    if not args.from_cache:
        run_compute(
            cancer=args.cancer,
            labels_path=args.labels,
            embedding_path=args.embedding,
            cache_dir=cache_dir,
            seed=args.seed,
            tsne_perplexity=args.tsne_perplexity,
        )

    written: List[str] = []
    for n in names:
        path = plotters[n](cache_dir=cache_dir, out_dir=out_dir, dpi=args.dpi)
        written.append(str(path))
        print("Wrote", path)

    from kirc_manuscript_cache_io import load_meta  # noqa: PLC0415

    meta = load_meta(cache_dir) if (cache_dir / "meta.json").is_file() else {}
    _write_readme(out_dir, written, meta)
    print("Done. Cache:", cache_dir)


if __name__ == "__main__":
    main()
