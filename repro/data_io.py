"""
Load multi-omics matrices the same way as Subtype-DCC (concatenate blocks, MinMax per paper pipeline).
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Iterable, Optional, Sequence, Set, Tuple

import numpy as np
import pandas as pd
from sklearn.preprocessing import MinMaxScaler

BLOCK_ORDER: Sequence[Tuple[str, str]] = (
    ("CN", "CN.fea"),
    ("meth", "meth.fea"),
    ("miRNA", "miRNA.fea"),
    ("rna", "rna.fea"),
)


def resolve_fea_root(repo_root: Optional[Path] = None) -> Path:
    """Feature matrices: prefer vendor/Subtype-DCC/subtype_file/fea (junction target)."""
    if repo_root is None:
        repo_root = Path(__file__).resolve().parents[1]
    vendored = repo_root / "vendor" / "Subtype-DCC" / "subtype_file" / "fea"
    if vendored.is_dir():
        return vendored
    return repo_root / "subtype_file" / "fea"


def patient_barcode(sample_id: str) -> str:
    s = str(sample_id).strip().upper()
    if len(s) >= 12 and s.startswith("TCGA"):
        return s[:12]
    return s


def load_concat_matrix(
    cancer_type: str,
    fea_root: Path,
    exclude_omics: Optional[Iterable[str]] = None,
    minmax: bool = True,
) -> Tuple[np.ndarray, list[str]]:
    """Returns X (n_samples, n_features), sample_ids as in matrix rows."""
    ex: Set[str] = set(exclude_omics or [])
    parts = []
    sample_names: Optional[list[str]] = None
    base = fea_root / cancer_type
    for key, fname in BLOCK_ORDER:
        if key in ex:
            continue
        path = base / fname
        if not path.is_file():
            raise FileNotFoundError(f"Missing feature file: {path}")
        df = pd.read_csv(path, header=0, index_col=0, sep=",")
        if sample_names is None:
            sample_names = list(df.columns)
        else:
            cols = list(df.columns)
            if cols != sample_names:
                missing = set(sample_names) - set(cols)
                extra = set(cols) - set(sample_names)
                if missing or extra:
                    raise ValueError(
                        f"{key} ({fname}): sample columns differ from first block "
                        f"(missing={len(missing)}, extra={len(extra)})."
                    )
                df = df.reindex(columns=sample_names)
        parts.append(df.values)
    if not parts:
        raise ValueError("No omics loaded.")
    X = np.concatenate(parts, axis=0).T.astype(np.float64)
    if minmax:
        X = MinMaxScaler().fit_transform(X)
    assert sample_names is not None
    return X, sample_names


def _canonical_header(name: str) -> str:
    """Match GerkeLab / cBioPortal headers: ignore case, collapse underscores vs spaces."""
    x = str(name).strip().lower().replace("_", " ")
    return re.sub(r"\s+", " ", x)


def _pick_col(df: pd.DataFrame, *candidates: str) -> Optional[pd.Series]:
    canon_map = {_canonical_header(c): c for c in df.columns}
    for cand in candidates:
        key = _canonical_header(cand)
        if key in canon_map:
            return df[canon_map[key]]
    return None


def _col_vals_aligned(n: int, s: Optional[pd.Series]) -> np.ndarray:
    """
    Build a length-n array from a column Series taken from ``sub``.
    Must use row order (``.values``), not index alignment: ``sub`` keeps the
    original dataframe integer index, while the output table uses
    ``patient_barcode`` as index — reindexing would otherwise yield all NaN.
    """
    if s is None:
        return np.array([np.nan] * n, dtype=object)
    v = s.to_numpy()
    if len(v) != n:
        raise ValueError(f"Clinical column length {len(v)} != patient rows {n}")
    return v


def _load_clinical_legacy(df: pd.DataFrame, cancer_type: str) -> pd.DataFrame:
    cols = [
        "bcr_patient_barcode",
        "acronym",
        "vital_status",
        "days_to_death",
        "days_to_last_followup",
        "gender",
        "age_at_initial_pathologic_diagnosis",
        "pathologic_M",
        "pathologic_N",
        "pathologic_T",
        "pathologic_stage",
    ]
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(f"Clinical file missing columns: {missing}")
    sub = df[cols].loc[df["acronym"] == cancer_type].copy()
    sub["patient_barcode"] = sub["bcr_patient_barcode"].astype(str).str.upper().str.slice(0, 12)
    sub["status"] = np.where(sub["vital_status"] == "Dead", 1, 0)
    sub["days"] = np.where(
        sub["status"] == 1,
        pd.to_numeric(sub["days_to_death"], errors="coerce"),
        pd.to_numeric(sub["days_to_last_followup"], errors="coerce"),
    )
    return sub.set_index("patient_barcode")


def _load_clinical_merged(df: pd.DataFrame, cancer_type: str) -> pd.DataFrame:
    """
    GerkeLab cBioportal_data.tsv / TCGA-CDR style wide tables (OS, OS.time, type, ...).
    Column names are matched case-insensitively where possible.
    """
    ct = cancer_type.upper()
    acronym_s = _pick_col(
        df,
        "tcga pancanatlas cancer type acronym",
        "acronym",
        "type",
        "cancer type abbreviation",
        "oncotree code",
    )
    if acronym_s is None:
        raise ValueError(
            "Could not find cancer-type column (tried acronym, type, ...). "
            "Use legacy clinical_PANCAN_patient_with_followup.tsv or check file headers."
        )
    sub = df.loc[acronym_s.astype(str).str.upper().str.strip() == ct].copy()
    if sub.empty:
        study = _pick_col(df, "study_id", "study id", "cancer study")
        if study is not None:
            pat = ct.lower()
            sub = df.loc[study.astype(str).str.lower().str.contains(pat, na=False)].copy()
    if sub.empty:
        raise ValueError(f"No rows for cancer type {cancer_type!r} after merge-table detection.")

    pid = _pick_col(
        sub,
        "patient id",
        "bcr patient barcode",
        "bcr_patient_barcode",
        "patient_id",
        "#patient identifier",
        "patient identifier",
    )
    if pid is None:
        raise ValueError("Could not find patient ID column (bcr_patient_barcode / PATIENT_ID).")
    sub["patient_barcode"] = pid.astype(str).str.upper().str.slice(0, 12)

    # Survival: prefer Liu OS endpoint (days); fallback vital_status + days columns
    os_time = _pick_col(
        sub,
        "overall survival (months)",
        "os.time",
        "os_time",
        "survival_time",
        "days_to_last_follow_up",
    )
    os_event = _pick_col(sub, "overall survival status", "os", "liu_os")
    if os_time is not None and os_event is not None:
        t = pd.to_numeric(os_time, errors="coerce")
        # cBioPortal/Gerke months -> approximate days (for axes comparable to TCGA day fields)
        time_col = os_time.name if hasattr(os_time, "name") else ""
        if time_col and "month" in str(time_col).lower():
            sub["days"] = t * (365.25 / 12.0)
        else:
            sub["days"] = t
        ev = os_event
        if pd.api.types.is_numeric_dtype(ev):
            sub["status"] = pd.to_numeric(ev, errors="coerce").fillna(0).astype(int)
        else:
            evs = ev.astype(str).str.upper()
            sub["status"] = np.where(
                evs.str.contains("DECEASED") | evs.str.contains("DEAD") | (evs == "1"),
                1,
                0,
            )
    else:
        vital = _pick_col(sub, "vital_status", "liu_vital_status")
        ddeath = _pick_col(sub, "days_to_death")
        dfoll = _pick_col(sub, "days_to_last_followup", "days_to_last_follow_up")
        if vital is None or ddeath is None or dfoll is None:
            raise ValueError(
                "Merged clinical file missing survival columns "
                "(need OS+OS.time or vital_status+days_to_death+days_to_last_followup)."
            )
        sub["status"] = np.where(vital.astype(str).str.lower().isin(("dead", "1")), 1, 0)
        sub["days"] = np.where(
            sub["status"] == 1,
            pd.to_numeric(ddeath, errors="coerce"),
            pd.to_numeric(dfoll, errors="coerce"),
        )

    sub = sub.drop_duplicates(subset=["patient_barcode"], keep="first")

    g = _pick_col(sub, "gender", "sex")
    age = _pick_col(sub, "age_at_initial_pathologic_diagnosis", "age", "diagnosis age")
    pm = _pick_col(
        sub,
        "american joint committee on cancer metastasis stage code",
        "pathologic_m",
        "pathologic_m_pathologic_spread",
    )
    pn = _pick_col(
        sub,
        "neoplasm disease lymph node stage american joint committee on cancer code",
        "pathologic_n",
        "pathologic_n_pathologic_spread",
    )
    pt = _pick_col(
        sub,
        "american joint committee on cancer tumor stage code",
        "pathologic_t",
        "pathologic_t_pathologic_spread",
    )
    ps = _pick_col(
        sub,
        "neoplasm disease stage american joint committee on cancer code",
        "pathologic_stage",
        "ajcc_pathologic_tumor_stage",
        "stage",
    )

    n = len(sub)
    age_arr = pd.to_numeric(_col_vals_aligned(n, age), errors="coerce")

    out = pd.DataFrame(
        {
            "gender": _col_vals_aligned(n, g),
            "age_at_initial_pathologic_diagnosis": age_arr,
            "pathologic_M": _col_vals_aligned(n, pm),
            "pathologic_N": _col_vals_aligned(n, pn),
            "pathologic_T": _col_vals_aligned(n, pt),
            "pathologic_stage": _col_vals_aligned(n, ps),
            "days": sub["days"].values,
            "status": sub["status"].values,
        },
        index=sub["patient_barcode"].values,
    )
    out.index.name = "patient_barcode"
    return out


def load_clinical_table(
    clinical_path: Path,
    cancer_type: str,
) -> pd.DataFrame:
    """
    Supports:
    - Legacy ``clinical_PANCAN_patient_with_followup.tsv`` (Broad/GDC supplemental).
    - GerkeLab ``cBioportal_data.tsv`` (download via scripts/download_clinical.py --source gerke).
    """
    df = pd.read_csv(clinical_path, sep="\t", header=0, low_memory=False)
    df.columns = [str(c).strip() for c in df.columns]

    legacy_need = {
        "bcr_patient_barcode",
        "acronym",
        "vital_status",
        "days_to_death",
        "days_to_last_followup",
    }
    if legacy_need.issubset(df.columns):
        return _load_clinical_legacy(df, cancer_type)
    return _load_clinical_merged(df, cancer_type)
