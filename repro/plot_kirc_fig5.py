"""Plot Fig 5 (GO BP/CC/MF bars) from cache only — no Enrichr calls."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPRO_DIR = Path(__file__).resolve().parent
if str(REPRO_DIR) not in sys.path:
    sys.path.insert(0, str(REPRO_DIR))

from kirc_manuscript_cache_io import default_cache_dir, load_go_fig5_table  # noqa: E402
from kirc_manuscript_lib import REPO_ROOT, fig5_go_bars_from_table  # noqa: E402


def plot(cache_dir: Path | None = None, out_dir: Path | None = None, dpi: int = 200) -> Path:
    cache_dir = cache_dir or default_cache_dir()
    out_dir = out_dir or (REPO_ROOT / "repro_out/kirc_manuscript_figs")
    go_df, uniq = load_go_fig5_table(cache_dir)
    if not uniq:
        raise FileNotFoundError(
            f"No GO table or cluster list (expected go_fig5.csv with cluster column). Cache: {cache_dir}"
        )
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / "KIRC_fig5_GO_bars.png"
    fig5_go_bars_from_table(go_df, uniq, out, dpi=dpi)
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
