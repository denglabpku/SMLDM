import numpy as np
import matplotlib.pyplot as plt
import torch
import torch.nn as nn
import torch.optim as optim
from torchvision import transforms
from torch.utils.data import DataLoader
import scipy.io as sio
import time
from pathlib import Path
# from skimage import io
import h5py
# from CNN_Model import buildModel, project_01, normalize_im
from ptorch_CNN_model import Model as buildModel, project_01, normalize_im
import argparse

def test_model(datafile, weights_file, meanstd_file, savename, upsampling_factor=8, debug=0):
    """
    This function tests a trained model on the desired test set, given the 
    tiff stack of test images, learned weights, and normalization factors.
    
    # Inputs
    datafile          - the tiff stack of test images 
    weights_file      - the saved weights file generated in train_model
    meanstd_file      - the saved mean and standard deviation file generated in train_model
    savename          - the filename for saving the recovered SR image
    upsampling_factor - the upsampling factor for reconstruction (default 8)
    debug             - boolean whether to save individual frame predictions (default 0)
    
    # Outputs
    function saves a mat file with the recovered image, and optionally saves 
    individual frame predictions in case debug=1. (default is debug=0)    
    """

    # # Load the tiff data
    # Images = io.imread(datafile)

    with h5py.File(datafile, 'r') as matfile:
        matfile = h5py.File(datafile, 'r')
        # Access the data
        cell_PSF = matfile['cell_PSF']
        xyI_dataset = cell_PSF[cell_PSF['xyI'][0][67]]  # Get specific dataset
        BoundxyI_dataset = cell_PSF[cell_PSF['BoundxyI'][0][67]]
        referenced_data = matfile[xyI_dataset[0, 1]]  # Obtain the data
        xyI_array = np.array(referenced_data)  # Convert to numpy array
        referenced_data = matfile[BoundxyI_dataset[0, 1]]  # Obtain the data
        BoundxyI_array = np.array(referenced_data)  # Convert to numpy array



        # Extract x, y, and intensity values
        x_values = xyI_array[0]
        y_values = xyI_array[1]
        intensity_values = xyI_array[2]

        # Find the brightest pixel
        brightest_pixel_index = np.argmax(intensity_values)

        # Create a 32x32 array
        new_array = np.zeros((32, 32))

        # Place the original array into the center of the new array
        x_offset = int(32 / 2) - int(x_values[brightest_pixel_index])
        y_offset = int(32 / 2) - int(y_values[brightest_pixel_index])

        for i in range(len(x_values)):
            new_array[int(x_values[i]) + x_offset, int(y_values[i]) + y_offset] = intensity_values[i]

        # Fill empty pixels with defined Gaussian noise (mean=0, std=1)
        mask = new_array == 0
        # noise = np.random.normal(0, 1, new_array.shape)
        noise = BoundxyI_array[2].mean()
        new_array[mask] = noise[mask]

        # # Optionally, apply gaussian filter to the array
        # new_array = gaussian_filter(new_array, sigma=1)

        # Now new_array contains the 32x32 array with the specified modifications
    

    # Get dataset dimensions
    K, M, N = Images.shape

    # Upsampling using a simple nearest neighbor interp.
    Images_upsampled = np.kron(Images, np.ones((1, upsampling_factor, upsampling_factor)))
    Images = Images_upsampled

    # Upsampled frames dimensions
    K, M, N = Images.shape

    # Build the model for a bigger image
    model = buildModel((M, N, 1))

    # Load the trained weights
    model.load_state_dict(torch.load(weights_file))

    # # Load mean and std
    # matfile = sio.loadmat(meanstd_file)
    # test_mean = matfile['mean_test']
    # test_std = matfile['std_test']

    # Setting type
    Images = Images.astype('float32')

    # Normalize each sample by its own mean and std
    Images_norm = np.zeros(Images.shape, dtype=np.float32)
    for i in range(Images.shape[0]):
        Images_norm[i, :, :] = project_01(Images[i, :, :])
        Images_norm[i, :, :] = normalize_im(Images_norm[i, :, :], test_mean, test_std)

    # Reshaping
    Images_norm = np.expand_dims(Images_norm, axis=3)

    # Make a prediction and time it
    start = time.time()
    predicted_density = model.predict(Images_norm, batch_size=1)
    end = time.time()
    print(end - start)

    # Threshold negative values
    predicted_density[predicted_density < 0] = 0

    # Resulting sum images
    WideField = np.squeeze(np.sum(Images_norm, axis=0))
    Recovery = np.squeeze(np.sum(predicted_density, axis=0))

    # Look at the sum image
    f, (ax1, ax2) = plt.subplots(1, 2, sharey=True, sharex=True)
    ax1.imshow(WideField)
    ax1.set_title('Wide Field')
    ax2.imshow(Recovery)
    ax2.set_title('Sum of Predictions')
    f.subplots_adjust(hspace=0)
    plt.show()

    # Save predictions to a matfile to open later in matlab
    mdict = {"Recovery": Recovery}
    sio.savemat(savename, mdict)

    # Save predicted density in each frame for debugging purposes
    if debug:
        mdict = {"Predictions": predicted_density}
        sio.savemat(savename + '_predictions.mat', mdict)

    return f


if __name__ == '__main__':

    # Start a parser
    parser = argparse.ArgumentParser()

    # Path of the tiff stack to be reconstructed
    parser.add_argument('--datafile', help="path to cell_PSF mat file")

    # Path of the optimal model weights and normalization factors, saved after training
    parser.add_argument('--weights_name', help="path to the trained model weights as pth-file")
    # parser.add_argument('--meanstd_name', help="path to the saved normalization factors as m-file")

    # Mean intensity of background
    parser.add_argument('--upsampling_factor', type=int, default=8, help="desired upsampling factor")

    # Standard deviation of background
    parser.add_argument('--upsampling_factor', type=int, default=8, help="desired upsampling factor")

    # Path for saving the Superresolution reconstruction matfile
    parser.add_argument('--savename', type=str, help="path for saving the Superresolution reconstruction matfile")

    # Upsampling factor for the superresolution reconstruction
    parser.add_argument('--upsampling_factor', type=int, default=8, help="desired upsampling factor")

    # Boolean debugging constant
    parser.add_argument('--debug', type=int, default=0, help="boolean (0/1) for saving individual predictions")

    # Parse the input arguments
    args = parser.parse_args()

    datafile = Path(args.datafile)
    weights_file = Path(args.weights_name)
    
    savename = Path(args.savename)
    upsampling_factor = args.upsampling_factor
    debug = args.debug

    # Run the testing/reconstruction process
    test_model(datafile, 
               weights_file,
               savename,
               upsampling_factor, 
               debug)
