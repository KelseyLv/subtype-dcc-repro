"""
Early-integration baselines on the same MinMax-scaled concatenated matrix as Subtype-DCC:
K-means, Spectral clustering, NMF + argmax cluster assignment (rank = n_clusters).
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans, SpectralClustering
from sklearn.decomposition import NMF

from data_io import load_concat_matrix, resolve_fea_root


def run_kmeans(X: np.ndarray, sample_names: list[str], k: int, seed: int, out: Path):
    model = KMeans(n_clusters=k, random_state=seed, n_init=10)
    lab = model.fit_predict(X) + 1
    pd.DataFrame({"sample_name": sample_names, "cluster": lab}).to_csv(out, sep="\t", index=False)


def run_spectral(X: np.ndarray, sample_names: list[str], k: int, seed: int, out: Path):
    model = SpectralClustering(n_clusters=k, random_state=seed, affinity="rbf", eigen_tol=1e-4)
    lab = model.fit_predict(X) + 1
    pd.DataFrame({"sample_name": sample_names, "cluster": lab}).to_csv(out, sep="\t", index=False)


def run_nmf(X: np.ndarray, sample_names: list[str], k: int, seed: int, out: Path):
    Xp = X - X.min() + 1e-9
    model = NMF(n_components=k, init="nndsvda", random_state=seed, max_iter=500)
    W = model.fit_transform(Xp)
    lab = np.argmax(W, axis=1) + 1
    pd.DataFrame({"sample_name": sample_names, "cluster": lab}).to_csv(out, sep="\t", index=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cancer", "-c", type=str, required=True)
    ap.add_argument("--method", type=str, choices=("kmeans", "spectral", "nmf"), required=True)
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--seed", type=int, default=21)
    ap.add_argument("--fea-root", type=Path, default=None)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    fea_root = args.fea_root or resolve_fea_root()
    X, samples = load_concat_matrix(args.cancer, fea_root, minmax=True)
    args.out.parent.mkdir(parents=True, exist_ok=True)

    if args.method == "kmeans":
        run_kmeans(X, samples, args.k, args.seed, args.out)
    elif args.method == "spectral":
        run_spectral(X, samples, args.k, args.seed, args.out)
    else:
        run_nmf(X, samples, args.k, args.seed, args.out)
    print("Wrote", args.out)


if __name__ == "__main__":
    main()
