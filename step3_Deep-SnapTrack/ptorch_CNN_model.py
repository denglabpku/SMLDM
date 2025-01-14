# Model definition + helper functions

import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
import numpy as np

# Define a function that projects an image to the range [0,1]
def project_01(im):
    im = np.squeeze(im)
    min_val = im.min()
    max_val = im.max()
    return (im - min_val) / (max_val - min_val)

# Normalize image given mean and std
def normalize_im(im, dmean, dstd):
    im = np.squeeze(im)
    im_norm = (im - dmean) / dstd
    return im_norm

# Define a MATLAB-like 2D Gaussian filter
def matlab_style_gauss2D(shape=(7,7), sigma=1):
    """ 
    2D gaussian filter - should give the same result as:
    MATLAB's fspecial('gaussian',[shape],[sigma]) 
    """
    m,n = [(ss-1.)/2. for ss in shape]
    y,x = np.ogrid[-m:m+1,-n:n+1]
    h = np.exp(-(x*x + y*y) / (2.*sigma*sigma))
    h[h < np.finfo(h.dtype).eps*h.max()] = 0
    sumh = h.sum()
    if sumh != 0:
        h /= sumh
    h = h*2.0
    h = h.astype(np.float32)
    return h

# Expanded the filter dimensions
psf_heatmap = matlab_style_gauss2D(shape=(7,7), sigma=1)
gfilter = torch.from_numpy(psf_heatmap).view(1, 1, 7, 7)

# Combined MSE + L1 loss
class L1L2loss(nn.Module):
    def __init__(self, input_shape, device):
        super(L1L2loss, self).__init__()
        self.input_shape = input_shape
        # Define gfilter as a module parameter
        self.gfilter = nn.Parameter(gfilter.to(device=device, dtype=torch.float32), requires_grad=False)        


    def forward(self, heatmap_true, spikes_pred):
        # Cast heatmap_true to float type to match the type of heatmap_pred
        heatmap_true = heatmap_true.float()

        # Cast spikes_pred to float32
        spikes_pred = spikes_pred.float()                

        # Generate the heatmap corresponding to the predicted spikes
        padding = (self.gfilter.size(-1) - 1) // 2  # Calculate padding size for 'same' padding
        spikes_pred_padded = F.pad(spikes_pred, (padding, padding, padding, padding), mode='constant', value=0)
        heatmap_pred = F.conv2d(spikes_pred_padded, self.gfilter)

        # Heatmaps MSE
        loss_heatmaps = F.mse_loss(heatmap_true, heatmap_pred)

        # L1 on the predicted spikes        
        loss_spikes = F.l1_loss(spikes_pred, torch.zeros_like(spikes_pred))
        
        return loss_heatmaps + loss_spikes
    
# Define the concatenated conv2, batch normalization, and relu block
    # rk and ck means row and column of kernel size
def conv_bn_relu(in_channels, nb_filter, rk, ck): 
    return nn.Sequential(
        nn.Conv2d(in_channels, nb_filter, kernel_size=(rk, ck), stride=1, padding="same", bias=False),
        nn.BatchNorm2d(nb_filter),
        nn.ReLU()
    )

# Define the model architecture
class CNN(nn.Module):
    def __init__(self):
        super(CNN, self).__init__()
        self.features1 = conv_bn_relu(1, 32, 3, 3)
        self.pool1 = nn.MaxPool2d(kernel_size=2, stride=2)
        self.features2 = conv_bn_relu(32, 64, 3, 3)
        self.pool2 = nn.MaxPool2d(kernel_size=2, stride=2)
        self.features3 = conv_bn_relu(64, 128, 3, 3)
        self.pool3 = nn.MaxPool2d(kernel_size=2, stride=2)
        self.features4 = conv_bn_relu(128, 512, 3, 3)
        self.up5 = nn.Upsample(scale_factor=2, mode='bilinear', align_corners=True)
        self.features5 = conv_bn_relu(512, 128, 3, 3)
        self.up6 = nn.Upsample(scale_factor=2, mode='bilinear', align_corners=True)
        self.features6 = conv_bn_relu(128, 64, 3, 3)
        self.up7 = nn.Upsample(scale_factor=2, mode='bilinear', align_corners=True)
        self.features7 = conv_bn_relu(64, 32, 3, 3)

    def forward(self, x):
        out = self.features1(x)
        out = self.pool1(out)
        out = self.features2(out)
        out = self.pool2(out)
        out = self.features3(out)
        out = self.pool3(out)
        out = self.features4(out)
        out = self.up5(out)
        out = self.features5(out)
        out = self.up6(out)
        out = self.features6(out)
        out = self.up7(out)
        out = self.features7(out)
        return out

# Define the model building for an arbitrary input size
class Model(nn.Module):
    def __init__(self):
        super(Model, self).__init__()
        self.cnn = CNN()
        self.density_pred = nn.Conv2d(32, 1, kernel_size=1, stride=1, padding=0)

    def forward(self, x):
        out = self.cnn(x)
        out = self.density_pred(out)
        return out





