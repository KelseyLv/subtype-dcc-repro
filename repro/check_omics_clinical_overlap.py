"""
Quick check: omics column IDs (from rna.fea header) vs clinical patient_barcode
for a cohort — same join logic as evaluate.py (no training required).

Usage (repo root):
  python repro/check_omics_clinical_overlap.py -c LIHC
  python repro/check_omics_clinical_overlap.py -c LIHC --clinical data/clinical/foo.tsv
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from data_io import load_clinical_table, patient_barcode, resolve_fea_root


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("-c", "--cancer", required=True, help="Cohort folder under fea/, e.g. LIHC")
    ap.add_argument(
        "--clinical",
        type=Path,
        default=None,
        help="Default: <repo>/data/clinical/cBioportal_data.tsv",
    )
    ap.add_argument("--rna-fea", type=Path, default=None, help="Override path to rna.fea")
    args = ap.parse_args()

    repo = Path(__file__).resolve().parents[1]
    clin_path = args.clinical or (repo / "data" / "clinical" / "cBioportal_data.tsv")
    if not clin_path.is_file():
        raise SystemExit(f"Clinical file not found: {clin_path}")

    fea = args.rna_fea or (resolve_fea_root(repo) / args.cancer.upper() / "rna.fea")
    if not fea.is_file():
        raise SystemExit(f"rna.fea not found: {fea}")

    df0 = pd.read_csv(fea, nrows=0, index_col=0)
    omics_cols = [str(x).strip() for x in df0.columns]
    omics_patients = {patient_barcode(s) for s in omics_cols}

    clin = load_clinical_table(clin_path, args.cancer.upper())
    clinical_patients = set(clin.index.astype(str))

    inter = omics_patients & clinical_patients
    only_o = omics_patients - clinical_patients
    only_c = clinical_patients - omics_patients

    lab = pd.DataFrame({"sample_name": omics_cols})
    lab["patient_barcode"] = lab["sample_name"].map(patient_barcode)
    merged = lab.join(clin, on="patient_barcode", how="inner")
    with_surv = merged.dropna(subset=["days", "status"])
    with_surv = with_surv[with_surv["days"] > 0]

    print(f"Cohort: {args.cancer.upper()}")
    print(f"rna.fea: {fea}")
    print(f"clinical: {clin_path}")
    print(f"Omics columns (aliquots/samples): {len(omics_cols)}")
    print(f"Unique patient_barcode from omics: {len(omics_patients)}")
    print(f"Clinical rows (patients): {len(clinical_patients)}")
    print(f"Patient intersection (omics and clinical): {len(inter)}")
    print(f"  Omics-only patients (no clinical row): {len(only_o)}")
    print(f"  Clinical-only patients (not in omics): {len(only_c)}")
    print(f"Simulated evaluate merge rows (per aliquot): {len(merged)}")
    print(f"Rows with valid survival (days, status non-null, days>0): {len(with_surv)}")
    # After train, evaluate also needs ≥2 clusters among merged rows; ID overlap is what we check here.
    ok_surv = len(with_surv) >= 2
    print(f"Enough rows for survival / clinical merge check (>=2 with survival): {ok_surv}")
    if only_o and len(only_o) <= 10:
        print(f"Omics-only examples: {sorted(only_o)}")
    elif only_o:
        print(f"Omics-only examples (5): {sorted(list(only_o))[:5]}")


if __name__ == "__main__":
    main()
