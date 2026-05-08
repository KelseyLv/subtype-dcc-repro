import torch
import numpy as np
import pandas as pd
from sklearn.preprocessing import MinMaxScaler
from torch.utils.data import DataLoader

BLOCK_FILES = [
    ("CN", "CN.fea"),
    ("meth", "meth.fea"),
    ("miRNA", "miRNA.fea"),
    ("rna", "rna.fea"),
]


def _fea_root():
    return "../subtype_file/fea/"


def load_omics_blocks(cancer_type, exclude_omics=None):
    """
    Load per-omics matrices (features x samples), stack features (same as original
    np.concatenate on axis=0 then .T -> samples x features).
    exclude_omics: set of keys among CN, meth, miRNA, rna
    """
    exclude_omics = exclude_omics or set()
    parts = []
    root = _fea_root() + cancer_type + "/"
    sample_names = None
    for key, fname in BLOCK_FILES:
        if key in exclude_omics:
            continue
        path = root + fname
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
        raise ValueError("All omics excluded; nothing to load.")
    feature = np.concatenate(parts, axis=0).T
    return feature, sample_names


def get_input_dim(cancer_type, exclude_omics=None):
    ex = set(exclude_omics) if exclude_omics else set()
    feat, _ = load_omics_blocks(cancer_type, exclude_omics=ex)
    return feat.shape[1]


def load_scaled_features(cancer_type, exclude_omics=None):
    """
    Full stacked matrix MinMax-scaled (same transform as training), as a single
    CPU float tensor (samples x features) plus sample column order from first block.
    """
    exclude_set = set(exclude_omics) if exclude_omics else set()
    feature, sample_names = load_omics_blocks(cancer_type, exclude_omics=exclude_set)
    scaler = MinMaxScaler()
    feature = scaler.fit_transform(feature)
    tensor = torch.tensor(feature, dtype=torch.float32)
    return tensor, sample_names


def get_feature(cancer_type, batch_size, training, exclude_omics=None, num_workers=0):
    exclude_set = set(exclude_omics) if exclude_omics else set()
    feature_tensor, _ = load_scaled_features(cancer_type, exclude_set)
    return DataLoader(
        feature_tensor,
        batch_size=batch_size,
        shuffle=training,
        num_workers=int(num_workers),
    )
