"""
Read/write the KIRC manuscript figure cache (heavy compute: t-SNE, Enrichr, Fig7 matrix).
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np
import pandas as pd

from kirc_manuscript_lib import CACHE_VERSION, REPO_ROOT


@dataclass
class KircManuscriptCache:
    meta: Dict[str, Any]
    Z: np.ndarray
    labels: np.ndarray
    samples: List[str]
    xlim: Tuple[float, float]
    ylim: Tuple[float, float]
    X_rna: np.ndarray
    genes: List[str]
    marker_idx: List[int]
    uniq_clusters: List[int]
    go_fig5: pd.DataFrame
    kegg_fig6: pd.DataFrame
    fig7_matrix: np.ndarray
    fig7_pathways: List[str]
    fig7_mirlabels: List[str]


def default_cache_dir() -> Path:
    return REPO_ROOT / "repro_out" / "kirc_manuscript_cache"


def _paths(cache_dir: Path) -> Dict[str, Path]:
    return {
        "meta": cache_dir / "meta.json",
        "Z": cache_dir / "Z.npy",
        "labels": cache_dir / "labels.npy",
        "X_rna": cache_dir / "X_rna.npy",
        "genes": cache_dir / "genes.json",
        "X_mir": cache_dir / "X_mir.npy",
        "mir_names": cache_dir / "mir_names.json",
        "go_fig5": cache_dir / "go_fig5.csv",
        "kegg_fig6": cache_dir / "kegg_fig6.csv",
        "fig7_matrix": cache_dir / "fig7_matrix.npy",
        "fig7_pathways": cache_dir / "fig7_pathways.json",
        "fig7_mirlabels": cache_dir / "fig7_mirlabels.json",
    }


def save_cache_bundle(
    cache_dir: Path,
    meta: Dict[str, Any],
    Z: np.ndarray,
    labels: np.ndarray,
    samples: List[str],
    xlim: Tuple[float, float],
    ylim: Tuple[float, float],
    X_rna: np.ndarray,
    genes: List[str],
    marker_idx: List[int],
    uniq_clusters: List[int],
    go_fig5: pd.DataFrame,
    kegg_fig6: pd.DataFrame,
    fig7_matrix: np.ndarray,
    fig7_pathways: List[str],
    fig7_mirlabels: List[str],
    X_mir: np.ndarray,
    mir_names: List[str],
) -> None:
    cache_dir.mkdir(parents=True, exist_ok=True)
    p = _paths(cache_dir)
    meta_out = {
        **meta,
        "cache_version": CACHE_VERSION,
        "samples": samples,
        "uniq_clusters": uniq_clusters,
        "marker_idx": marker_idx,
        "xlim": list(xlim),
        "ylim": list(ylim),
    }
    p["meta"].write_text(json.dumps(meta_out, indent=2), encoding="utf-8")
    np.save(p["Z"], Z.astype(np.float64))
    np.save(p["labels"], labels.astype(np.int32))
    np.save(p["X_rna"], X_rna.astype(np.float64))
    p["genes"].write_text(json.dumps(genes), encoding="utf-8")
    np.save(p["X_mir"], X_mir.astype(np.float64))
    p["mir_names"].write_text(json.dumps(mir_names), encoding="utf-8")
    go_fig5.to_csv(p["go_fig5"], index=False)
    kegg_fig6.to_csv(p["kegg_fig6"], index=False)
    np.save(p["fig7_matrix"], fig7_matrix.astype(np.float64))
    p["fig7_pathways"].write_text(json.dumps(fig7_pathways), encoding="utf-8")
    p["fig7_mirlabels"].write_text(json.dumps(fig7_mirlabels), encoding="utf-8")


def load_cache(cache_dir: Path) -> KircManuscriptCache:
    p = _paths(cache_dir)
    if not p["meta"].is_file():
        raise FileNotFoundError(f"Missing cache meta: {p['meta']}")
    meta = json.loads(p["meta"].read_text(encoding="utf-8"))
    Z = np.load(p["Z"])
    labels = np.load(p["labels"])
    samples = list(meta["samples"])
    xlim = tuple(float(x) for x in meta["xlim"])
    ylim = tuple(float(y) for y in meta["ylim"])
    X_rna = np.load(p["X_rna"])
    genes = json.loads(p["genes"].read_text(encoding="utf-8"))
    marker_idx = [int(x) for x in meta["marker_idx"]]
    uniq_clusters = [int(x) for x in meta["uniq_clusters"]]
    go_fig5 = pd.read_csv(p["go_fig5"]) if p["go_fig5"].is_file() else pd.DataFrame()
    kegg_fig6 = pd.read_csv(p["kegg_fig6"]) if p["kegg_fig6"].is_file() else pd.DataFrame()
    fig7_matrix = np.load(p["fig7_matrix"]) if p["fig7_matrix"].is_file() else np.array([[np.nan]])
    fig7_pathways = json.loads(p["fig7_pathways"].read_text(encoding="utf-8")) if p["fig7_pathways"].is_file() else []
    fig7_mirlabels = json.loads(p["fig7_mirlabels"].read_text(encoding="utf-8")) if p["fig7_mirlabels"].is_file() else []
    return KircManuscriptCache(
        meta=meta,
        Z=Z,
        labels=labels,
        samples=samples,
        xlim=xlim,
        ylim=ylim,
        X_rna=X_rna,
        genes=genes,
        marker_idx=marker_idx,
        uniq_clusters=uniq_clusters,
        go_fig5=go_fig5,
        kegg_fig6=kegg_fig6,
        fig7_matrix=fig7_matrix,
        fig7_pathways=fig7_pathways,
        fig7_mirlabels=fig7_mirlabels,
    )


def cache_is_complete(cache_dir: Path) -> bool:
    p = _paths(cache_dir)
    required = ["meta", "Z", "labels", "X_rna", "genes", "go_fig5", "kegg_fig6", "fig7_matrix"]
    return all(p[k].is_file() for k in required)


def load_meta(cache_dir: Path) -> Dict[str, Any]:
    return json.loads((cache_dir / "meta.json").read_text(encoding="utf-8"))


def load_tsne_bundle(cache_dir: Path) -> Tuple[np.ndarray, np.ndarray, Tuple[float, float], Tuple[float, float]]:
    meta = load_meta(cache_dir)
    Z = np.load(_paths(cache_dir)["Z"])
    labels = np.load(_paths(cache_dir)["labels"])
    xlim = tuple(float(x) for x in meta["xlim"])
    ylim = tuple(float(y) for y in meta["ylim"])
    return Z, labels, xlim, ylim


def load_rna_marker_bundle(
    cache_dir: Path,
) -> Tuple[np.ndarray, List[str], List[int], List[int], np.ndarray]:
    meta = load_meta(cache_dir)
    X_rna = np.load(_paths(cache_dir)["X_rna"])
    genes: List[str] = json.loads(_paths(cache_dir)["genes"].read_text(encoding="utf-8"))
    marker_idx = [int(x) for x in meta["marker_idx"]]
    uniq_clusters = [int(x) for x in meta["uniq_clusters"]]
    labels = np.load(_paths(cache_dir)["labels"])
    return X_rna, genes, marker_idx, uniq_clusters, labels


def load_go_fig5_table(cache_dir: Path) -> Tuple[pd.DataFrame, List[int]]:
    p = _paths(cache_dir)["go_fig5"]
    if not p.is_file():
        meta = load_meta(cache_dir)
        return pd.DataFrame(), [int(x) for x in meta.get("uniq_clusters", [])]
    df = pd.read_csv(p)
    if df.empty or "cluster" not in df.columns:
        meta = load_meta(cache_dir)
        return df, [int(x) for x in meta.get("uniq_clusters", [])]
    df["cluster"] = df["cluster"].astype(int)
    uniq = sorted(int(x) for x in df["cluster"].unique())
    return df, uniq


def load_kegg_fig6_table(cache_dir: Path) -> pd.DataFrame:
    p = _paths(cache_dir)["kegg_fig6"]
    return pd.read_csv(p) if p.is_file() else pd.DataFrame()


def load_fig7_bundle(cache_dir: Path) -> Tuple[np.ndarray, List[str], List[str]]:
    p = _paths(cache_dir)
    mat = np.load(p["fig7_matrix"]) if p["fig7_matrix"].is_file() else np.array([[np.nan]])
    pathways = json.loads(p["fig7_pathways"].read_text(encoding="utf-8")) if p["fig7_pathways"].is_file() else []
    mirlabels = json.loads(p["fig7_mirlabels"].read_text(encoding="utf-8")) if p["fig7_mirlabels"].is_file() else []
    return mat, pathways, mirlabels
