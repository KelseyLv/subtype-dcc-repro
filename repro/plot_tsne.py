"""t-SNE on embedding .fea from Subtype-DCC (Figure 4A style)."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.manifold import TSNE

from evaluate import load_labels_dcc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fea", type=Path, required=True, help="Subtype-DCC results/{CANCER}.fea")
    ap.add_argument("--labels", type=Path, required=True, help=".dcc file")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--perplexity", type=float, default=30.0)
    ap.add_argument("--seed", type=int, default=21)
    args = ap.parse_args()

    emb = pd.read_csv(args.fea, sep="\t", index_col=0)
    lab = load_labels_dcc(args.labels)
    lab = lab.set_index("sample_name")
    common = emb.index.intersection(lab.index)
    X = emb.loc[common].values.astype(np.float32)
    y = lab.loc[common, "cluster"].values

    tsne = TSNE(
        n_components=2,
        perplexity=min(args.perplexity, max(5, X.shape[0] // 4)),
        random_state=args.seed,
        init="random",
        learning_rate=max(50.0, X.shape[0] / 12.0),
    )
    Z = tsne.fit_transform(X)

    plt.figure(figsize=(7, 6))
    for c in np.unique(y):
        m = y == c
        plt.scatter(Z[m, 0], Z[m, 1], s=12, label=f"C{int(c)}", alpha=0.7)
    plt.legend(title="Subtype", markerscale=2)
    plt.xlabel("t-SNE 1")
    plt.ylabel("t-SNE 2")
    plt.title("t-SNE of Subtype-DCC embedding")
    plt.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(args.out, dpi=150)
    plt.close()
    print("Saved", args.out)


if __name__ == "__main__":
    main()
