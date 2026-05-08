import os
import sys
import argparse
from pathlib import Path

import numpy as np
import torch
from modules import ae, network, contrastive_loss
from utils import yaml_config_hook, save_checkpoint_last, save_model
from torch.utils.data import DataLoader

from dataloader import get_input_dim, load_scaled_features
import pandas as pd
from tensorboardX import SummaryWriter
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def inference(loader, model, device):
    model.eval()
    cluster_vector = []
    feature_vector = []
    for _, x in enumerate(loader):
        x = x.float().to(device)
        with torch.no_grad():
            c, h = model.forward_cluster(x)
            cluster_vector.extend(c.cpu().detach().numpy())
            feature_vector.extend(h.cpu().detach().numpy())
    cluster_vector = np.array(cluster_vector)
    feature_vector = np.array(feature_vector)
    print("Features shape {}".format(feature_vector.shape))
    return cluster_vector, feature_vector


def train_one_epoch(model, DL, optimizer, device, loss_device, cluster_number, args):
    loss_epoch = 0
    for _, x in enumerate(DL):
        optimizer.zero_grad()
        x_i = (x + torch.normal(0, 1, size=(x.shape[0], x.shape[1]))).float().to(device)
        x_j = (x + torch.normal(0, 1, size=(x.shape[0], x.shape[1]))).float().to(device)
        z_i, z_j, c_i, c_j = model(x_i, x_j)
        criterion_instance = contrastive_loss.DCL(temperature=0.5, weight_fn=None)
        criterion_cluster = contrastive_loss.ClusterLoss(
            cluster_number, args.cluster_temperature, loss_device
        ).to(loss_device)
        loss_instance = criterion_instance(z_i, z_j) + criterion_instance(z_j, z_i)
        loss_cluster = criterion_cluster(c_i, c_j)
        loss = loss_instance + loss_cluster
        loss.backward()
        optimizer.step()
        loss_epoch += loss.item()
    return loss_epoch


def _repo_root() -> Path:
    """train.py lives in <repo>/vendor/Subtype-DCC/Subtype-DCC/."""
    return Path(__file__).resolve().parents[3]


def _ensure_repro_on_path() -> None:
    repro = _repo_root() / "repro"
    s = str(repro)
    if s not in sys.path:
        sys.path.insert(0, s)


def _resolve_clinical_path(cli_value: str) -> Path | None:
    """Absolute path to clinical TSV; None if missing."""
    repo = _repo_root()
    if cli_value and str(cli_value).strip():
        p = Path(cli_value)
        if not p.is_absolute():
            p = repo / p
        return p if p.is_file() else None
    default = repo / "data" / "clinical" / "cBioportal_data.tsv"
    return default if default.is_file() else None


def _argmax_cluster_labels(model, loader, device) -> np.ndarray:
    model.eval()
    parts = []
    with torch.no_grad():
        for x in loader:
            x = x.float().to(device)
            c, _ = model.forward_cluster(x)
            parts.append(c.cpu().numpy().astype(np.int64, copy=False))
    return np.concatenate(parts) if parts else np.array([], dtype=np.int64)


def _neg_log10_logrank_from_model(
    model,
    loader,
    device,
    sample_names: list,
    clinical_df: pd.DataFrame,
) -> tuple[float, int]:
    """
    Same survival contrast as repro/evaluate.py (multivariate log-rank).
    Returns (neg_log10_p, n_merged); nan neg_log10 if not computable.
    """
    _ensure_repro_on_path()
    from evaluate import merge_labels_clinical, survival_neg_log10_p

    y = _argmax_cluster_labels(model, loader, device)
    if y.size == 0 or len(sample_names) != len(y):
        return float("nan"), 0
    lab = pd.DataFrame(
        {"sample_name": list(sample_names), "cluster": (y + 1).astype(int)},
    )
    merged = merge_labels_clinical(lab, clinical_df, cluster_col="cluster")
    if merged.empty or merged["cluster"].nunique() < 2:
        return float("nan"), int(len(merged))
    return float(survival_neg_log10_p(merged)), int(len(merged))


def _cluster_assignment_stats(model, loader, device, cluster_number):
    """Argmax cluster counts on full data (loader must be shuffle=False, same scaling as train)."""
    model.eval()
    labels = []
    with torch.no_grad():
        for x in loader:
            x = x.float().to(device)
            c, _ = model.forward_cluster(x)
            labels.append(c.cpu().numpy().astype(np.int64, copy=False))
    y = np.concatenate(labels) if labels else np.array([], dtype=np.int64)
    counts = np.bincount(y, minlength=cluster_number)
    total = int(counts.sum())
    p = counts.astype(np.float64) / max(total, 1)
    nonempty = int((counts > 0).sum())
    ent = float(-np.sum(np.where(p > 0, p * np.log(p + 1e-12), 0.0)))
    h_max = float(np.log(cluster_number)) if cluster_number > 1 else 1.0
    return counts, nonempty, ent, h_max


def _loss_window_rel_range(window):
    """(max-min)/|mean| on last W losses; small => flat / stable."""
    if len(window) < 2:
        return 1.0
    a = np.asarray(window, dtype=np.float64)
    m = float(np.mean(a))
    if abs(m) < 1e-12:
        return 0.0
    return float((np.max(a) - np.min(a)) / (abs(m) + 1e-12))


def draw_fig(losses, name):
    x1 = range(0, len(losses))
    save_file = "./results/" + name + "Train_loss.png"
    os.makedirs(os.path.dirname(save_file), exist_ok=True)
    plt.cla()
    plt.title("Train loss vs. epoch", fontsize=20)
    plt.plot(x1, losses, ".-")
    plt.xlabel("epoch", fontsize=20)
    plt.ylabel("Train loss", fontsize=20)
    plt.grid()
    plt.savefig(save_file)
    plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    parser.add_argument("--cancer_type", "-c", type=str, default="BRCA")
    parser.add_argument("--batch_size", type=int, default=64)
    parser.add_argument("--cluster_number", type=int, default=-1)
    parser.add_argument(
        "--exclude_omics",
        type=str,
        default="",
        help="Comma-separated omics to drop: CN,meth,miRNA,rna (ablation)",
    )

    cancer_dict = {
        "BRCA": 5,
        "BLCA": 5,
        "KIRC": 4,
        "LUAD": 3,
        "PAAD": 2,
        "SKCM": 4,
        "STAD": 3,
        "UCEC": 4,
        "UVM": 4,
        "GBM": 2,
    }

    config = yaml_config_hook("config/config.yaml")
    for k, v in config.items():
        parser.add_argument(f"--{k}", default=v, type=type(v))
    parser.add_argument(
        "--log_every",
        type=int,
        default=10,
        help="Print training loss every N epochs (default: 10).",
    )
    parser.add_argument(
        "--cluster_log_every",
        type=int,
        default=0,
        help="If >0, every N epochs print argmax cluster counts + balance entropy (0=off).",
    )
    parser.add_argument(
        "--clinical",
        type=str,
        default="",
        help="Clinical TSV for mid-training log-rank (default: <repo>/data/clinical/cBioportal_data.tsv).",
    )
    parser.add_argument(
        "--survival_logrank_every",
        type=int,
        default=10,
        help="Every N epochs compute survival neg_log10(log-rank p) like evaluate.py (0=off).",
    )
    parser.add_argument(
        "--no_early_stop",
        action="store_true",
        help="Disable plateau / patience early stopping; always run for --epochs.",
    )
    parser.add_argument(
        "--plateau_window",
        type=int,
        default=30,
        help="Number of past epochs to measure loss flatness (max-min)/mean.",
    )
    parser.add_argument(
        "--plateau_rel_tol",
        type=float,
        default=0.005,
        help="If (max-min)/mean in the window is below this, treat as 'flat' (default 0.005 = 0.5%%).",
    )
    parser.add_argument(
        "--plateau_stable_epochs",
        type=int,
        default=40,
        help="Stop after this many consecutive 'flat' epochs (default 40).",
    )
    parser.add_argument(
        "--min_epochs",
        type=int,
        default=100,
        help="No early stop before this many epochs (default 100).",
    )
    parser.add_argument(
        "--improvement_patience",
        type=int,
        default=0,
        help="If >0, stop when best loss does not improve by --improvement_min_delta for this many epochs (0=off).",
    )
    parser.add_argument(
        "--improvement_min_delta",
        type=float,
        default=1.0,
        help="Minimal absolute loss drop to count as improvement for improvement_patience.",
    )
    args = parser.parse_args()

    cluster_number = (
        cancer_dict[args.cancer_type]
        if args.cluster_number < 0
        else args.cluster_number
    )
    print("Clusters:", cluster_number)

    exclude_set = {s.strip() for s in args.exclude_omics.split(",") if s.strip()}

    model_path = "./save/" + args.cancer_type
    os.makedirs(model_path, exist_ok=True)
    os.makedirs("./results", exist_ok=True)

    torch.manual_seed(args.seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(args.seed)
        torch.cuda.manual_seed(args.seed)
    np.random.seed(args.seed)

    logger = SummaryWriter(log_dir="./log")

    input_dim = get_input_dim(args.cancer_type, exclude_omics=exclude_set)
    full_tensor, sample_name_list = load_scaled_features(
        args.cancer_type, exclude_omics=exclude_set
    )
    nw = int(getattr(args, "workers", 0) or 0)
    DL = DataLoader(
        full_tensor,
        batch_size=args.batch_size,
        shuffle=True,
        num_workers=nw,
    )
    cluster_eval_loader = DataLoader(
        full_tensor,
        batch_size=args.batch_size,
        shuffle=False,
        num_workers=nw,
    )

    clinical_df = None
    survival_every = int(getattr(args, "survival_logrank_every", 0) or 0)
    if survival_every > 0:
        cpath = _resolve_clinical_path(getattr(args, "clinical", "") or "")
        if cpath is None:
            print(
                "[log-rank] disabled: pass --clinical or place data/clinical/cBioportal_data.tsv under repo.",
                flush=True,
            )
            survival_every = 0
        else:
            _ensure_repro_on_path()
            try:
                from data_io import load_clinical_table

                clinical_df = load_clinical_table(cpath, args.cancer_type)
                print(
                    f"[log-rank] enabled every {survival_every} epochs using {cpath}",
                    flush=True,
                )
            except Exception as e:
                print(f"[log-rank] disabled (clinical load failed): {e}", flush=True)
                survival_every = 0

    ae_module = ae.AE(input_dim=input_dim)
    model = network.Network(ae_module, args.feature_dim, cluster_number)
    model = model.to(device)

    optimizer = torch.optim.Adam(
        model.parameters(), lr=args.learning_rate, weight_decay=args.weight_decay
    )
    loss_device = device

    losses = []
    early_stop = not args.no_early_stop
    plateau_stable_run = 0
    best_loss = float("inf")
    imp_wait = 0
    finished_epoch = args.start_epoch - 1

    for epoch in range(args.start_epoch, args.epochs + 1):
        model.train()
        loss_epoch = train_one_epoch(
            model, DL, optimizer, device, loss_device, cluster_number, args
        )
        losses.append(loss_epoch)
        finished_epoch = epoch
        logger.add_scalar("train loss", loss_epoch, epoch)

        if epoch % args.log_every == 0 or epoch == args.epochs:
            rel = _loss_window_rel_range(losses[-args.plateau_window :])
            print(
                f"Epoch [{epoch}/{args.epochs}]\t Loss: {loss_epoch:.6f}\t "
                f"window_rel_range({min(len(losses), args.plateau_window)}): {rel:.6f}",
                flush=True,
            )

        if args.cluster_log_every > 0 and (
            epoch % args.cluster_log_every == 0 or epoch == args.epochs
        ):
            counts, nz, ent, h_max = _cluster_assignment_stats(
                model, cluster_eval_loader, device, cluster_number
            )
            ratio = ent / h_max if h_max > 1e-12 else 0.0
            print(
                f"  [clusters] epoch {epoch}: counts={counts.tolist()} "
                f"nonempty={nz}/{cluster_number} H/Hmax={ratio:.4f}",
                flush=True,
            )
            logger.add_scalar("cluster_nonempty", nz, epoch)
            logger.add_scalar("cluster_balance_H_ratio", ratio, epoch)
            for k in range(cluster_number):
                logger.add_scalar(f"cluster_count/k{k}", int(counts[k]), epoch)

        if survival_every > 0 and clinical_df is not None and (
            epoch % survival_every == 0 or epoch == args.epochs
        ):
            try:
                nlp, n_ev = _neg_log10_logrank_from_model(
                    model,
                    cluster_eval_loader,
                    device,
                    sample_name_list,
                    clinical_df,
                )
                if np.isfinite(nlp):
                    print(
                        f"  [log-rank] epoch {epoch}: neg_log10_p={nlp:.6f} "
                        f"n_eval={n_ev} (multivariate log-rank, same as evaluate.py)",
                        flush=True,
                    )
                    logger.add_scalar("survival_neg_log10_logrank_p", nlp, epoch)
                else:
                    print(
                        f"  [log-rank] epoch {epoch}: skipped (n_eval={n_ev}, need >=2 clusters with overlap).",
                        flush=True,
                    )
            except Exception as e:
                print(f"  [log-rank] epoch {epoch}: failed: {e}", flush=True)

        save_checkpoint_last(model_path, model, optimizer, epoch)
        save_model(model_path, model, optimizer, args.epochs)

        if early_stop and epoch >= args.min_epochs:
            if len(losses) >= args.plateau_window:
                w = losses[-args.plateau_window :]
                rel_w = _loss_window_rel_range(w)
                if rel_w < args.plateau_rel_tol:
                    plateau_stable_run += 1
                else:
                    plateau_stable_run = 0
                if plateau_stable_run >= args.plateau_stable_epochs:
                    print(
                        f"\n[Early stop] Loss plateau: relative range < {args.plateau_rel_tol} "
                        f"for {args.plateau_stable_epochs} consecutive epochs (epoch {epoch}).\n",
                        flush=True,
                    )
                    break

            if args.improvement_patience > 0:
                if loss_epoch < best_loss - args.improvement_min_delta:
                    best_loss = loss_epoch
                    imp_wait = 0
                else:
                    imp_wait += 1
                if imp_wait >= args.improvement_patience:
                    print(
                        f"\n[Early stop] No improvement >= {args.improvement_min_delta} "
                        f"for {args.improvement_patience} epochs (epoch {epoch}).\n",
                        flush=True,
                    )
                    break

    draw_fig(losses, args.cancer_type)
    print(
        f"Training finished at epoch {finished_epoch}/{args.epochs} "
        f"(final loss sum over batches: {losses[-1]:.6f}).",
        flush=True,
    )

    dataloader = DataLoader(
        full_tensor,
        batch_size=args.batch_size,
        shuffle=False,
        num_workers=nw,
    )

    model_eval = network.Network(
        ae.AE(input_dim=input_dim),
        args.feature_dim,
        cluster_number,
    )
    model_fp = os.path.join(model_path, "checkpoint_last.tar")
    ckpt = torch.load(model_fp, map_location=device.type)
    print(f"Loaded checkpoint: {model_fp} (saved epoch {ckpt.get('epoch', '?')})", flush=True)
    model_eval.load_state_dict(ckpt["net"])
    model_eval.to(device)

    print("### Creating features from model ###")
    X, h = inference(dataloader, model_eval, device)
    output = pd.DataFrame(columns=["sample_name", "dcc"])
    # Row order of `X` / `h` matches `load_omics_blocks` (first block in BLOCK_FILES = CN),
    # not necessarily `rna.fea` on-disk column order.
    sample_name = list(sample_name_list)
    output["sample_name"] = sample_name
    output["dcc"] = X + 1
    out_file = "./results/" + args.cancer_type + ".dcc"
    output.to_csv(out_file, index=False, sep="\t")

    fea_out_file = "./results/" + args.cancer_type + ".fea"
    fea = pd.DataFrame(
        data=h,
        index=sample_name,
        columns=map(lambda x: "v" + str(x), range(h.shape[1])),
    )
    fea.to_csv(fea_out_file, header=True, index=True, sep="\t")
