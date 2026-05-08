import os
import torch


def save_model(model_path, model, optimizer, current_epoch):
    out = os.path.join(model_path, "checkpoint_{}.tar".format(current_epoch))
    state = {
        "net": model.state_dict(),
        "optimizer": optimizer.state_dict(),
        "epoch": current_epoch,
    }
    torch.save(state, out)


def save_checkpoint_last(model_path, model, optimizer, epoch):
    """Latest weights for inference (supports early stopping)."""
    out = os.path.join(model_path, "checkpoint_last.tar")
    torch.save(
        {
            "net": model.state_dict(),
            "optimizer": optimizer.state_dict(),
            "epoch": epoch,
        },
        out,
    )
