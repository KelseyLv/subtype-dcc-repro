"""
Launch Subtype-DCC training with one omics layer removed (paper Figure S1 style).
Run from repo root: python repro/run_ablation.py --cancer KIRC
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

OMIC_KEYS = ("CN", "meth", "miRNA", "rna")
TRAIN_DIR = Path(__file__).resolve().parents[1] / "vendor" / "Subtype-DCC" / "Subtype-DCC"
TRAIN_PY = TRAIN_DIR / "train.py"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cancer", "-c", type=str, required=True)
    ap.add_argument(
        "--drop",
        type=str,
        choices=OMIC_KEYS,
        action="append",
        default=[],
        help="Omic to drop (repeat flag for multiple). Default: run each single-omic removal.",
    )
    ap.add_argument("--extra-args", type=str, default="", help="Appended to train.py (e.g. --epochs 50)")
    args = ap.parse_args()

    to_run = args.drop if args.drop else list(OMIC_KEYS)
    for key in to_run:
        cmd = [
            sys.executable,
            str(TRAIN_PY),
            "-c",
            args.cancer,
            "--exclude_omics",
            key,
        ]
        if args.extra_args:
            cmd.extend(args.extra_args.split())
        print("Running:", " ".join(cmd))
        subprocess.run(cmd, cwd=TRAIN_DIR, check=True)


if __name__ == "__main__":
    main()
