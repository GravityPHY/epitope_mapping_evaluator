import argparse
import os
import random
import torch
import torch.nn.parameter
import torch.nn.parallel
import torch.optim as optim
import torch.utils.data
from models.PointNet.model import PointNetDenseCls12, feature_transform_regularizer
from models.PointNet.dataset import ShapeNetDataset3aug

import torch.multiprocessing
torch.multiprocessing.set_sharing_strategy('file_system')

def pairwise_distances(x, y=None):
    '''
    Input: x is a Nxd matrix
           y is an optional Mxd matirx
    Output: dist is a NxM matrix where dist[i,j] is the square norm between x[i,:] and y[j,:]
            if y is not given then use 'y=x'.
    i.e. dist[i,j] = ||x[i,:]-y[j,:]||^2
    '''
    x_norm = (x ** 2).sum(1).view(-1, 1)
    if y is not None:
        y_t = torch.transpose(y, 0, 1)
        y_norm = (y ** 2).sum(1).view(1, -1)
    else:
        y_t = torch.transpose(x, 0, 1)
        y_norm = x_norm.view(1, -1)

    dist = x_norm + y_norm - 2.0 * torch.mm(x, y_t)
    return torch.clamp(dist, 0.0, np.inf)

def gk(x):
    cen=torch.nn.Parameter(torch.tensor([0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95]).float(), requires_grad=False).cuda()
    xmat=x.float().cuda()-cen
    sigma=200
    y = torch.sum(torch.sigmoid(sigma * (xmat + 0.1 / 2)) - torch.sigmoid(sigma * (xmat - 0.1 / 2)),dim=0)
    y = y / torch.sum(y)
    return y

