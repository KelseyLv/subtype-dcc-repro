"""
Paper metrics: -log10(log-rank p), count of significant clinical variables (6 tests).
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Dict, Tuple

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, kruskal

from data_io import load_clinical_table, patient_barcode

try:
    from lifelines.statistics import multivariate_logrank_test
except ImportError:
    multivariate_logrank_test = None


CLINICAL_COLS = [
    ("age_at_initial_pathologic_diagnosis", "numeric"),
    ("gender", "categorical"),
    ("pathologic_stage", "categorical"),
    ("pathologic_T", "categorical"),
    ("pathologic_N", "categorical"),
    ("pathologic_M", "categorical"),
]


def merge_labels_clinical(
    labels: pd.DataFrame,
    clinical: pd.DataFrame,
    sample_col: str = "sample_name",
    cluster_col: str = "cluster",
) -> pd.DataFrame:
    labels = labels.copy()
    labels["patient_barcode"] = labels[sample_col].map(patient_barcode)
    merged = labels.join(clinical, on="patient_barcode", how="inner")
    return merged


def survival_neg_log10_p(merged: pd.DataFrame, cluster_col: str = "cluster") -> float:
    if multivariate_logrank_test is None:
        raise RuntimeError("Install lifelines: pip install lifelines")
    df = merged.dropna(subset=["days", "status", cluster_col]).copy()
    df = df[df["days"] > 0]
    if df[cluster_col].nunique() < 2:
        return 0.0
    r = multivariate_logrank_test(
        event_durations=df["days"],
        groups=df[cluster_col],
        event_observed=df["status"],
    )
    p = float(r.p_value)
    if p <= 0:
        return 300.0
    return -math.log10(p)


def clinical_significant_count(
    merged: pd.DataFrame,
    cluster_col: str = "cluster",
    alpha: float = 0.05,
) -> Tuple[int, Dict[str, float]]:
    pvals: Dict[str, float] = {}
    sig = 0
    for col, kind in CLINICAL_COLS:
        sub = merged[[cluster_col, col]].dropna()
        if sub.empty or sub[cluster_col].nunique() < 2:
            pvals[col] = 1.0
            continue
        if kind == "numeric":
            groups = [g[col].astype(float).values for _, g in sub.groupby(cluster_col)]
            groups = [g[~np.isnan(g)] for g in groups if len(g) > 0]
            if len(groups) < 2:
                pvals[col] = 1.0
                continue
            _, p = kruskal(*groups)
        else:
            ct = pd.crosstab(sub[cluster_col], sub[col])
            if ct.size == 0 or ct.shape[0] < 2 or ct.shape[1] < 2:
                pvals[col] = 1.0
                continue
            _, p, _, _ = chi2_contingency(ct)
        p = float(p)
        pvals[col] = p
        if p < alpha:
            sig += 1
    return sig, pvals


def load_labels_dcc(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    if "dcc" in df.columns:
        df = df.rename(columns={"dcc": "cluster"})
    return df


def main():
    ap = argparse.ArgumentParser(description="Evaluate cluster labels vs survival & clinical (paper metrics).")
    ap.add_argument("--labels", type=Path, required=True, help="TSV with sample_name, cluster (or dcc)")
    ap.add_argument("--clinical", type=Path, required=True, help="clinical_PANCAN_patient_with_followup.tsv")
    ap.add_argument("--cancer", "-c", type=str, required=True)
    ap.add_argument("--out-json", type=Path, default=None)
    args = ap.parse_args()

    clin = load_clinical_table(args.clinical, args.cancer)
    labels = load_labels_dcc(args.labels)
    merged = merge_labels_clinical(labels, clin)
    if merged.empty:
        raise SystemExit("No overlap between labels and clinical barcodes.")

    nlp = survival_neg_log10_p(merged)
    n_sig, pvals = clinical_significant_count(merged)
    out = {
        "n_samples_eval": int(len(merged)),
        "neg_log10_logrank_p": nlp,
        "n_significant_clinical": n_sig,
        "clinical_pvals": pvals,
    }
    print(json.dumps(out, indent=2))
    if args.out_json:
        args.out_json.parent.mkdir(parents=True, exist_ok=True)
        args.out_json.write_text(json.dumps(out, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
