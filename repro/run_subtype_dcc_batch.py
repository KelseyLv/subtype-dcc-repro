"""Train Subtype-DCC for all nine TCGA cohorts (requires subtype_file/fea populated)."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

CANCERS_9 = ["BRCA", "BLCA", "KIRC", "LUAD", "PAAD", "SKCM", "STAD", "UCEC", "UVM"]
TRAIN_DIR = Path(__file__).resolve().parents[1] / "vendor" / "Subtype-DCC" / "Subtype-DCC"
TRAIN_PY = TRAIN_DIR / "train.py"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cancers", nargs="+", default=CANCERS_9)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--extra-args", type=str, default="", help="e.g. --epochs 10 for smoke test")
    args = ap.parse_args()

    for c in args.cancers:
        cmd = [sys.executable, str(TRAIN_PY), "-c", c]
        if args.extra_args:
            cmd.extend(args.extra_args.split())
        print("[Subtype-DCC]", " ".join(cmd))
        if args.dry_run:
            continue
        subprocess.run(cmd, cwd=TRAIN_DIR, check=True)


if __name__ == "__main__":
    main()
