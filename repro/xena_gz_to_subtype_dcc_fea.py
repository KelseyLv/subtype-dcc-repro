"""
Convert UCSC Xena `.gz` (tab-separated) downloads into Subtype-DCC `*.fea` CSVs.

Xena exports: rows = features, columns = samples (first column = feature id).
Subtype-DCC expects: same layout, but **comma**-separated, filenames
`rna.fea`, `miRNA.fea`, `meth.fea`, `CN.fea` under `subtype_file/fea/<CANCER>/`.

Methylation matrices are huge (~450k probes): this script uses **chunked reads**
+ a min-heap to pick top-variance probes without loading the full matrix into RAM.

Usage (repo root)::

    .venv\\Scripts\\python.exe repro\\xena_gz_to_subtype_dcc_fea.py ^
      --cancer LIHC ^
      --dir vendor\\Subtype-DCC\\subtype_file\\fea\\LIHC ^
      --rna rna-TCGA.LIHC.sampleMap_HiSeqV2.gz ^
      --mirna miRNA-TCGA.LIHC.sampleMap_miRNA_HiSeq_gene.gz ^
      --meth meth-TCGA.LIHC.sampleMap_HumanMethylation450.gz ^
      --cn CN-TCGA.LIHC.sampleMap_Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz
"""

from __future__ import annotations

import argparse
import heapq
from pathlib import Path
from typing import Dict, List, Set, Tuple

import numpy as np
import pandas as pd


def _header_sample_cols(path: Path) -> List[str]:
    """Sample aliquot IDs only (exclude first column = feature index)."""
    df0 = pd.read_csv(path, sep="\t", compression="gzip", nrows=0)
    return [str(c).strip() for c in df0.columns[1:]]


def _intersect_sample_cols(paths: Dict[str, Path]) -> List[str]:
    common: Set[str] | None = None
    for name, p in paths.items():
        cols = _header_sample_cols(p)
        s = set(cols)
        common = s if common is None else (common & s)
    if not common:
        raise ValueError("No common sample columns across all four matrices.")
    return sorted(common)


def _index_col_name(path: Path) -> str:
    df0 = pd.read_csv(path, sep="\t", compression="gzip", nrows=0)
    return str(df0.columns[0]).strip()


def _read_matrix_usecols(path: Path, usecols: List[str]) -> pd.DataFrame:
    df = pd.read_csv(
        path,
        sep="\t",
        index_col=0,
        header=0,
        compression="gzip",
        usecols=usecols,
        low_memory=False,
    )
    return df.apply(pd.to_numeric, errors="coerce")


def _fillna_row_median(df: pd.DataFrame) -> pd.DataFrame:
    df = df.replace([np.inf, -np.inf], np.nan)
    row_med = df.median(axis=1)
    return df.T.fillna(row_med).T


def _top_variance_rows(df: pd.DataFrame, max_rows: int) -> pd.DataFrame:
    if max_rows <= 0 or len(df) <= max_rows:
        return df
    v = np.nanvar(df.values.astype(np.float64), axis=1, ddof=0)
    order = np.argsort(-v)
    return df.iloc[order[:max_rows]]


def _meth_top_probe_ids(
    path: Path,
    sample_cols: List[str],
    k: int,
    chunksize: int = 25000,
) -> List[str]:
    """Two-phase: stream chunks, keep k probes with largest row-wise variance."""
    index_name = _index_col_name(path)
    cols_arg = [index_name] + sample_cols

    heap: List[Tuple[float, str]] = []
    for chunk in pd.read_csv(
        path,
        sep="\t",
        index_col=0,
        header=0,
        compression="gzip",
        usecols=cols_arg,
        chunksize=chunksize,
        low_memory=False,
    ):
        block = chunk.apply(pd.to_numeric, errors="coerce")
        vars_ = np.nanvar(block.values.astype(np.float64), axis=1, ddof=0)
        for probe, var in zip(block.index.astype(str), vars_):
            vv = float(var) if np.isfinite(var) else 0.0
            if len(heap) < k:
                heapq.heappush(heap, (vv, probe))
            elif vv > heap[0][0]:
                heapq.heapreplace(heap, (vv, probe))
    return [pid for _, pid in sorted(heap, key=lambda x: -x[0])]


def _meth_subset_from_ids(path: Path, sample_cols: List[str], probe_ids: List[str]) -> pd.DataFrame:
    want = set(probe_ids)
    index_name = _index_col_name(path)
    cols_arg = [index_name] + sample_cols
    parts: List[pd.DataFrame] = []
    for chunk in pd.read_csv(
        path,
        sep="\t",
        index_col=0,
        header=0,
        compression="gzip",
        usecols=cols_arg,
        chunksize=40000,
        low_memory=False,
    ):
        sub = chunk.loc[chunk.index.astype(str).isin(want)]
        if len(sub):
            parts.append(sub)
    if not parts:
        raise RuntimeError("Methylation subset empty.")
    out = pd.concat(parts, axis=0)
    out = out.loc[[p for p in probe_ids if p in out.index]]
    return out.apply(pd.to_numeric, errors="coerce")


def main() -> None:
    ap = argparse.ArgumentParser(description="Xena .gz -> Subtype-DCC *.fea (comma CSV).")
    ap.add_argument("--cancer", required=True, help="Cohort folder name, e.g. LIHC")
    ap.add_argument("--dir", type=Path, required=True, help="Directory containing the four .gz files")
    ap.add_argument("--rna", type=str, required=True, help="RNA matrix filename inside --dir")
    ap.add_argument("--mirna", type=str, required=True)
    ap.add_argument("--meth", type=str, required=True)
    ap.add_argument("--cn", type=str, required=True)
    ap.add_argument("--rna-max-genes", type=int, default=3500)
    ap.add_argument("--mirna-max", type=int, default=450)
    ap.add_argument("--meth-max-probes", type=int, default=3139)
    ap.add_argument("--cn-max-genes", type=int, default=3500)
    ap.add_argument(
        "--primary-tumor-only",
        action="store_true",
        help="Keep columns whose ID ends with -01A / -01B / -01 (heuristic primary).",
    )
    args = ap.parse_args()

    # Do not call `.resolve()` here: on Windows a junction under `vendor/.../fea` may
    # dereference to `external/...`; writing via the repo-relative path keeps intent clear.
    d = args.dir if args.dir.is_absolute() else (Path.cwd() / args.dir)
    paths = {
        "rna": d / args.rna,
        "mirna": d / args.mirna,
        "meth": d / args.meth,
        "cn": d / args.cn,
    }
    for k, p in paths.items():
        if not p.is_file():
            raise FileNotFoundError(f"Missing {k}: {p}")

    print("Scanning headers for common samples...")
    sample_order = _intersect_sample_cols(paths)
    if args.primary_tumor_only:
        sample_order = [
            c
            for c in sample_order
            if str(c).strip().upper().endswith(("01A", "01B", "01"))
        ]
    print(f"Using {len(sample_order)} samples (columns).")

    usecols_rna = [_index_col_name(paths["rna"])] + sample_order
    usecols_mir = [_index_col_name(paths["mirna"])] + sample_order
    usecols_cn = [_index_col_name(paths["cn"])] + sample_order
    print("Loading RNA / miRNA / CN (subset columns)...")
    rna = _read_matrix_usecols(paths["rna"], usecols_rna)
    mir = _read_matrix_usecols(paths["mirna"], usecols_mir)
    cn = _read_matrix_usecols(paths["cn"], usecols_cn)

    print("Pass 1/2 methylation: streaming variance to pick top probes...")
    top_meth_probes = _meth_top_probe_ids(paths["meth"], sample_order, args.meth_max_probes)
    print("Pass 2/2 methylation: loading selected probes...")
    meth = _meth_subset_from_ids(paths["meth"], sample_order, top_meth_probes)

    aligned = {
        "rna": _fillna_row_median(rna),
        "mirna": _fillna_row_median(mir),
        "meth": _fillna_row_median(meth),
        "cn": _fillna_row_median(cn),
    }

    rna_o = _top_variance_rows(aligned["rna"], args.rna_max_genes)
    mir_o = _top_variance_rows(aligned["mirna"], args.mirna_max)
    meth_o = aligned["meth"]
    cn_o = _top_variance_rows(aligned["cn"], args.cn_max_genes)

    out_map = [
        ("rna", rna_o, "rna.fea"),
        ("miRNA", mir_o, "miRNA.fea"),
        ("meth", meth_o, "meth.fea"),
        ("CN", cn_o, "CN.fea"),
    ]
    d.mkdir(parents=True, exist_ok=True)
    for _, mat, fname in out_map:
        dest = d / fname
        mat.to_csv(dest, index=True, header=True, sep=",")
        print(f"Wrote {dest} shape={mat.shape}")

    meta = d / "xena_to_fea_meta.txt"
    meta.write_text(
        f"cancer={args.cancer}\n"
        f"n_samples={len(sample_order)}\n"
        f"rna_rows={len(rna_o)} mirna_rows={len(mir_o)} meth_rows={len(meth_o)} cn_rows={len(cn_o)}\n"
        f"caps: rna_max={args.rna_max_genes} mirna_max={args.mirna_max} "
        f"meth_max={args.meth_max_probes} cn_max={args.cn_max_genes}\n",
        encoding="utf-8",
    )
    print("Wrote", meta)


if __name__ == "__main__":
    main()
