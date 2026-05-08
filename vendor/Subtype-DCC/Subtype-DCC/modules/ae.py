import torch.nn as nn


def block(in_c, out_c):
    layers = [nn.Linear(in_c, out_c), nn.ReLU(True)]
    return layers


class Encoder(nn.Module):
    def __init__(self, input_dim=9844, inter_dims=None):
        super(Encoder, self).__init__()
        if inter_dims is None:
            inter_dims = [5000, 2000, 1000, 256]
        self.encoder = nn.Sequential(
            nn.Dropout(),
            *block(input_dim, inter_dims[0]),
            *block(inter_dims[0], inter_dims[1]),
            *block(inter_dims[1], inter_dims[2]),
            *block(inter_dims[2], inter_dims[3]),
        )

    def forward(self, x):
        return self.encoder(x)


class Decoder(nn.Module):
    def __init__(self, inter_dims=None):
        super(Decoder, self).__init__()
        if inter_dims is None:
            inter_dims = [5000, 2000, 1000, 256]
        self.decoder = nn.Sequential(
            *block(inter_dims[-1], inter_dims[-2]),
            *block(inter_dims[-2], inter_dims[-3]),
            *block(inter_dims[-3], inter_dims[-4]),
        )

    def forward(self, z):
        return self.decoder(z)


class AE(nn.Module):
    def __init__(self, hid_dim=256, input_dim=9844, inter_dims=None):
        super(AE, self).__init__()
        if inter_dims is None:
            inter_dims = [5000, 2000, 1000, hid_dim]
        else:
            inter_dims = list(inter_dims)
            inter_dims[-1] = hid_dim
        self.encoder = Encoder(input_dim=input_dim, inter_dims=inter_dims)
        self.decoder = Decoder(inter_dims=inter_dims)
        self.rep_dim = hid_dim

    def forward(self, x):
        z = self.encoder(x)
        return z
