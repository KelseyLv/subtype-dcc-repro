"""Plot Fig 4A from `kirc_manuscript_compute_cache.py` output only (no Enrichr)."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPRO_DIR = Path(__file__).resolve().parent
if str(REPRO_DIR) not in sys.path:
    sys.path.insert(0, str(REPRO_DIR))

from kirc_manuscript_cache_io import default_cache_dir, load_tsne_bundle  # noqa: E402
from kirc_manuscript_lib import REPO_ROOT, fig4a_tsne  # noqa: E402


def plot(cache_dir: Path | None = None, out_dir: Path | None = None, dpi: int = 200) -> Path:
    cache_dir = cache_dir or default_cache_dir()
    out_dir = out_dir or (REPO_ROOT / "repro_out/kirc_manuscript_figs")
    Z, labels, xlim, ylim = load_tsne_bundle(cache_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / "KIRC_fig4A_tsne.png"
    fig4a_tsne(Z, labels, out, dpi=dpi, xlim=xlim, ylim=ylim)
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
