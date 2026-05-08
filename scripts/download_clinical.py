"""
Download Pan-Cancer clinical resources.

1) UCSC Xena Pan-Cancer Atlas phenotype (dense) — public HTTP, large TSV.
2) Optionally rename / symlink to clinical_PANCAN_patient_with_followup.tsv for Subtype-GAN scripts.

Our repro/evaluate.py expects columns like the legacy PANCAN TSV. Xena uses long format;
run `python repro/adapt_xena_clinical.py` after download to build the wide table, OR use
`--source gerke` to fetch GerkeLab merged file (smaller, already wide-style clinical variables).

Default: GerkeLab cBioportal_data.tsv (~9MB) — includes survival and many clinical fields.
"""

from __future__ import annotations

import argparse
import sys
import urllib.request
from pathlib import Path

GERKE_CBIO = (
    "https://raw.githubusercontent.com/GerkeLab/TCGAclinical/master/data/cBioportal_data.tsv"
)
XENA_PHENO = "https://pancanatlas.xenahubs.net/download/TCGA_phenotype_denseDataOnlyDownload.tsv"


def download(url: str, dest: Path, label: str) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"Downloading {label} ...")
    print(f"  URL: {url}")
    print(f"  -> {dest}")

    def reporthook(blocknum, blocksize, totalsize):
        if totalsize > 0 and blocknum % 128 == 0:
            got = blocknum * blocksize
            pct = min(100, 100 * got / totalsize)
            print(f"  ... {pct:.1f}%", end="\r", file=sys.stderr)

    urllib.request.urlretrieve(url, dest, reporthook=reporthook)
    print("\nDone.", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--source",
        choices=("gerke", "xena"),
        default="gerke",
        help="gerke: smaller merged TSV from GerkeLab/TCGAclinical; xena: raw phenotype matrix (very large)",
    )
    ap.add_argument(
        "--out",
        type=Path,
        default=Path("data/clinical/cBioportal_data.tsv"),
        help="Output file path",
    )
    args = ap.parse_args()

    if args.source == "gerke":
        download(GERKE_CBIO, args.out, "GerkeLab TCGAclinical")
    else:
        out = args.out
        if out.suffix != ".tsv":
            out = out.with_suffix(".tsv")
        download(XENA_PHENO, out, "Xena TCGA phenotype")


if __name__ == "__main__":
    main()
