"""
Snapshot metrics while Subtype-DCC is still training (no need to wait for .dcc).

Reads ``save/<CANCER>/checkpoint_last.tar`` under ``vendor/Subtype-DCC/Subtype-DCC``,
runs argmax clustering on the same MinMax-scaled matrix as ``train.py``, then prints:
  - checkpoint epoch
  - per-cluster sample counts (1-based labels, same convention as ``.dcc``)
  - ``neg_log10_logrank_p`` and ``n_significant_clinical`` (same logic as ``repro/evaluate.py``)

Run from **repository root** (so ``repro/data_io.py`` resolves). Training can keep running;
on Windows, avoid reading the checkpoint at the exact instant the trainer is writing it
(wait until a fresh ``Epoch [...]`` line just appeared).

Example (LIHC, match your train flags)::

  .venv\\Scripts\\python.exe repro\\monitor_subtype_dcc_checkpoint.py -c LIHC ^
    --cluster_number 4 --batch_size 16 --workers 0 ^
    --clinical data\\clinical\\cBioportal_data.tsv
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _sdcc_root() -> Path:
    return _repo_root() / "vendor" / "Subtype-DCC" / "Subtype-DCC"


def main() -> None:
    ap = argparse.ArgumentParser(description="Mid-training snapshot: clusters + log-rank metrics.")
    ap.add_argument("--cancer", "-c", type=str, required=True)
    ap.add_argument(
        "--checkpoint",
        type=str,
        default="",
        help="Path to checkpoint_last.tar (default: <Subtype-DCC>/save/<cancer>/checkpoint_last.tar).",
    )
    ap.add_argument(
        "--clinical",
        type=Path,
        default=Path("data/clinical/cBioportal_data.tsv"),
        help="Clinical TSV (same as evaluate.py).",
    )
    ap.add_argument("--cluster_number", type=int, required=True)
    ap.add_argument("--feature_dim", type=int, default=128, help="Must match config used in training.")
    ap.add_argument("--batch_size", type=int, default=64, help="Inference minibatch size only.")
    ap.add_argument("--workers", type=int, default=0, help="DataLoader workers for inference.")
    ap.add_argument(
        "--exclude_omics",
        type=str,
        default="",
        help="Comma-separated: CN,meth,miRNA,rna (must match training).",
    )
    ap.add_argument("--out-json", type=Path, default=None)
    args = ap.parse_args()

    repo = _repo_root()
    sdcc = _sdcc_root()
    if not sdcc.is_dir():
        raise SystemExit(f"Missing Subtype-DCC code dir: {sdcc}")

    ckpt_path = Path(args.checkpoint) if args.checkpoint else sdcc / "save" / args.cancer / "checkpoint_last.tar"
    if not ckpt_path.is_file():
        raise SystemExit(f"Checkpoint not found: {ckpt_path.resolve()}")

    clinical_path = args.clinical if args.clinical.is_absolute() else (repo / args.clinical).resolve()
    if not clinical_path.is_file():
        raise SystemExit(f"Clinical file not found: {clinical_path}")

    exclude_set = {s.strip() for s in args.exclude_omics.split(",") if s.strip()}

    repro_dir = str(repo / "repro")
    if repro_dir not in sys.path:
        sys.path.insert(0, repro_dir)
    from data_io import load_clinical_table
    from evaluate import clinical_significant_count, merge_labels_clinical, survival_neg_log10_p

    old_cwd = os.getcwd()
    os.chdir(sdcc)
    sys.path.insert(0, str(sdcc))
    try:
        import numpy as np
        import pandas as pd
        import torch
        from torch.utils.data import DataLoader

        from dataloader import get_input_dim, load_scaled_features
        from modules import ae, network

        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        input_dim = get_input_dim(args.cancer, exclude_omics=exclude_set)
        full_tensor, sample_names = load_scaled_features(args.cancer, exclude_omics=exclude_set)
        loader = DataLoader(
            full_tensor,
            batch_size=args.batch_size,
            shuffle=False,
            num_workers=int(args.workers),
        )

        model = network.Network(
            ae.AE(input_dim=input_dim),
            args.feature_dim,
            args.cluster_number,
        ).to(device)
        ckpt = torch.load(ckpt_path, map_location=device.type)
        model.load_state_dict(ckpt["net"], strict=True)
        ckpt_epoch = int(ckpt.get("epoch", -1))

        model.eval()
        labels = []
        with torch.no_grad():
            for x in loader:
                x = x.float().to(device)
                c, _ = model.forward_cluster(x)
                labels.append(c.cpu().numpy().astype(np.int64, copy=False))
        y = np.concatenate(labels)
        counts = np.bincount(y, minlength=args.cluster_number).tolist()
        labels_1 = (y + 1).astype(np.int64)
    finally:
        os.chdir(old_cwd)
        if sdcc in sys.path:
            sys.path.remove(sdcc)

    clin = load_clinical_table(clinical_path, args.cancer)
    lab_df = pd.DataFrame(
        {"sample_name": list(sample_names), "cluster": labels_1.astype(int)},
    )
    merged = merge_labels_clinical(lab_df, clin, cluster_col="cluster")
    if merged.empty:
        raise SystemExit("No overlap between inferred labels and clinical barcodes.")

    nlp = survival_neg_log10_p(merged)
    n_sig, pvals = clinical_significant_count(merged)
    if nlp <= 0:
        logrank_p = None
    elif nlp >= 299.0:
        logrank_p = 0.0
    else:
        logrank_p = float(10 ** (-nlp))

    out = {
        "cancer": args.cancer,
        "checkpoint": str(ckpt_path.resolve()),
        "checkpoint_epoch": ckpt_epoch,
        "cluster_argmax_counts_0based_order": counts,
        "n_samples_in_fea": int(len(sample_names)),
        "n_samples_eval": int(len(merged)),
        "logrank_p": logrank_p,
        "neg_log10_logrank_p": nlp,
        "n_significant_clinical": n_sig,
        "clinical_pvals": pvals,
    }
    print(json.dumps(out, indent=2, ensure_ascii=False))
    if args.out_json:
        args.out_json = args.out_json if args.out_json.is_absolute() else (repo / args.out_json)
        args.out_json.parent.mkdir(parents=True, exist_ok=True)
        args.out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
        print(f"Wrote {args.out_json}", flush=True)


if __name__ == "__main__":
    main()
