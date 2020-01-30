import torch.nn as nn
import sparseconvnet as scn
import torch.nn.functional as F

full_scale = 4096
residual_blocks=False
block_reps=1
m = 16
scale=80

class Model(nn.Module):
    def __init__(self, Feature_channel, num_class):
        nn.Module.__init__(self)
        dimension = 3
        self.compress_fn = nn.Linear(Feature_channel,dimension)
        self.bn = nn.BatchNorm1d(3)
        self.sparseModel = scn.Sequential().add(
           scn.InputLayer(dimension,full_scale, mode=4)).add(
           scn.SubmanifoldConvolution(dimension, 3, m, 3, False)).add(
               scn.UNet(dimension, block_reps, [m, 2*m, 3*m, 4*m, 5*m, 6*m, 7*m], residual_blocks)).add(
           scn.BatchNormReLU(m)).add(
           scn.OutputLayer(dimension))
        self.linear = nn.Linear(m, num_class)
    def forward(self,x):
        x[1] = self.bn(F.relu(self.compress_fn(x[1])))
        x=self.sparseModel(x)
        x=self.linear(x)
        return x

