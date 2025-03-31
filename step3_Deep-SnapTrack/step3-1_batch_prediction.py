import h5py
import nd2
import numpy as np
import torch
from os import path, listdir
import pandas as pd
from tqdm import tqdm
import torch.nn.functional as F
from ptorch_CNN_model import Model as buildModel
from skimage.measure import regionprops

def IndividualNormalize(img):    
    # Define a function that projects an image to the range [0,1]
    img = np.squeeze(img)
    min_val = img.min()
    max_val = img.max()
    img = (img - min_val) / (max_val - min_val)

    # Normalize image given mean and std
    mean_img = img.mean()
    std_img = img.std()
    img = (img - mean_img) / std_img    

    # img = img.reshape(1, img.shape[0], img.shape[1])  # Reshape to (1, height, width)

    return img

def simple_mask_cutoff(pred_img, min_threshold=0.1):
    # rescale image so that in the range [0 1]
    rescale_pred_img = (pred_img-pred_img.min())/(pred_img.max()-pred_img.min())

    # Find the index where the cumulative distribution first exceeds 
    # use contrast (Max and Min)
    pixel_threshold_index = rescale_pred_img >= min_threshold

    mask = np.zeros(pred_img.shape)

    # Get the corresponding pixel value
    mask[pixel_threshold_index] = 1

    area = mask.sum()

    return mask, area

def estimate_background_image(ND2_img,device):

    # device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    mean_values = np.mean(ND2_img, axis=(1, 2))
    normalized_img = ND2_img / mean_values[:, np.newaxis, np.newaxis]
    bug_frameIdx =  np.where(mean_values == 0)[0] # identify bug frame in ND2 (all-zero frame)
    normalized_img[bug_frameIdx,] = normalized_img[bug_frameIdx-1,] # use the next frame value, which should be valid frame
    normalized_img = torch.from_numpy(normalized_img).to(device)

    kernel_size = 401
    gap_size = 50
    normalized_img = F.pad(normalized_img.permute(1,2,0),(kernel_size // 2,kernel_size // 2),mode='reflect')

    normalized_img = normalized_img.permute(2,0,1)

    # Unfold the tensor to create overlapping windows
    normalized_img = normalized_img.unfold(0, kernel_size, gap_size)  # This will create a tensor of shape (1900, 256, 256, 101)

    # Reshape to apply median across the windows
    normalized_img = normalized_img.permute(3,0,1,2).reshape(kernel_size, -1)

    # Apply median along the first dimension
    median_filtered_values, _ = torch.median(normalized_img, dim=0)

    # Reshape the median filtered tensor back to its original shape
    output_tensor = median_filtered_values.view(-1, ND2_img.shape[1], ND2_img.shape[2])

    output_tensor = F.interpolate(output_tensor.permute(1,2,0), scale_factor=gap_size, mode='linear')
    output_tensor = output_tensor.permute(2,0,1).cpu().numpy()

    # Step 4: Rescale the median values by multiplying with the mean frame intensity
    bk_img = output_tensor[:mean_values.shape[0],:,:] * mean_values[:,np.newaxis,np.newaxis]

    return bk_img

def estimate_D(SRArea,coeff_a=1162.0):
    D = (SRArea/coeff_a)**2
    return D

# >>>>>>>>>>>>>>>>>>>>>>>>>> input >>>>>>>>>>>>>>>>>>>>>>>>>> #

weights_file = './checkpoints/MBX_20231220_110nmPix_rep2_epoch9.pth'

# parent directory where you save UNet motion blur extracted mat folder
rootDir = '/path/to/save/results/'

# directory of UNet motion blur extracted mat folder
dataDir = [   
    '20240712_Clust01_U2OS_Paxillin_Cell01',    
    ]

# Number of ND2 images per sample in dataDir
# ND2PerSample = [10 10 10 7 10 10 9]

UNet_model = 'UNet_mask_MBX_20240620_epoch20_Ch1'
blurmat_file_prefix = 'Blurdata_'+UNet_model
fitresult_file = 'Fitresult_'+UNet_model+'.csv'
sr_fileName = UNet_model+'_SR_pred_v3.csv'

# raw image files, keep order same as dataDir
ND2File = [    
    '/path/to/your/data/20240712_Clust01_U2OS_Paxillin_30p5ms_2kframe_001.nd2',    
    ]

# <<<<<<<<<<<<<<<<<<<<<<<<<<<< input <<<<<<<<<<<<<<<<<<<<<<<<<<<< #

# --------------- main code ------------ #
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(f'Using device {device}')

# Build the model
model = buildModel().to(device)

# Load the trained weights
model.load_state_dict(torch.load(weights_file))
model.eval()
print(f'Model loaded from {weights_file}')

ND2_index_adjusted = 0

for iSamp, nowDir in enumerate(dataDir):
    torch.cuda.empty_cache()
    # Create a dictionary with the column names as keys and empty lists as values
    data = {
        'Frame': [], 'PSF_idx': [], 'Xpos': [], 'Ypos': [], 'TotalPhoton': [], 'Intensity': [],
        'Background': [],'EllipticityIdx': [],'Angle': [], 'SNR': [],'COV': [],'UNetArea': [], 'SRArea': [], 'D': []
    }

    # Create a DataFrame from the dictionary
    df = pd.DataFrame(data)
    csv_save_path = path.join(rootDir,nowDir)
    df.to_csv(path.join(csv_save_path, sr_fileName), index=False, mode='w')   

    # datafile = path.join(rootDir,nowDir,blurmat_file)

    tab_psf_fitresult_path = path.join(rootDir,nowDir,fitresult_file)

    tab_psf_fitresult = pd.read_csv(tab_psf_fitresult_path)

    irow = 0
    iFrame = 0

    print(f'Processing file {nowDir}')

    files = listdir(path.join(rootDir, nowDir))
    sorted_files = sorted(files)  # Sort files alphabetically by name

    for file in sorted_files:
        if blurmat_file_prefix in file:            
            datafile = path.join(rootDir, nowDir, file)

            with h5py.File(datafile, 'r') as matfile:
                # matfile = h5py.File(datafile, 'r')
                # Access the data
                cell_PSF = matfile['cell_PSF']
                # total_frame_num = cell_PSF['xyI'][0].size
                # total_frame_num = int(matfile['impars']['tot_img_num'][0][0])

                # bk_img = []

                # for iSlice in range(ND2PerSample[iSamp]):

                ND2_img = nd2.imread(ND2File[ND2_index_adjusted])

                #----- version2: robust background estimation ---- #
                #----- see Marten Postma Sci Rep paper --# 
                slice_interval = 1000  
        
                # Calculate the number of full slices we can get  
                num_full_slices = ND2_img.shape[0] // slice_interval  
                        
                # Iterate over the full slices  
                for loop_idx, iSlice in enumerate(range(0, num_full_slices * slice_interval, slice_interval)):  
                    frame_slice = ND2_img[iSlice:iSlice + slice_interval, :, :]  
                    # Perform your analysis on frame_slice  
                    print(f"Shape of full frame slice: {frame_slice.shape}")  
                    bk_img = estimate_background_image(frame_slice,device)
                    # bk_img = temp_bk_img                
                    torch.cuda.empty_cache()

                    for isubFrame in tqdm(range(slice_interval)):
                    
                        xyI_dataset = cell_PSF[cell_PSF['xyI'][0][isubFrame+loop_idx*slice_interval]]  # last index is frame number
                            
                        total_PSF_num = xyI_dataset[0].size
                        if xyI_dataset.len() == 1: # this frame contains PSF, if 2 means empty frame       

                            for iPSF in range(total_PSF_num):
                            
                                referenced_data = matfile[xyI_dataset[0, iPSF]]  # last index is the PSF index in the current frame
                                xyI_array = np.array(referenced_data)  # Convert to numpy array
                                UNetArea = xyI_array.shape[1] # Store the U-Net area

                                # Extract x, y, and intensity values
                                x_values = xyI_array[0]
                                y_values = xyI_array[1]
                                intensity_values = xyI_array[2]
                                
                                center_column = (x_values.min()+x_values.max())/2
                                center_row = (y_values.min()+y_values.max())/2
                                    
                                # Place the original array into the center of the new array
                                x_offset = int(32 / 2) - round(center_column)
                                y_offset = int(32 / 2) - round(center_row)

                                # Create a 32x32 array
                                Image = np.zeros((32, 32))

                                # Check if out-of-Image dimension, discard if it is too big, most like UNet mask cover more than 2 molecules
                                if (y_values+y_offset).min() >= 0 and (y_values+y_offset).max() < 32 and (x_values+x_offset).min() >= 0 and (x_values+x_offset).max() < 32:

                                    # Get background and noise from bk_img
                                    background = np.mean(bk_img[isubFrame,y_values.astype(int)-1,x_values.astype(int)-1])
                                    noise = np.std(bk_img[isubFrame,y_values.astype(int)-1,x_values.astype(int)-1])
                                    
                                    for i in range(len(x_values)):
                                        Image[int(y_values[i]) + y_offset,int(x_values[i]) + x_offset] = intensity_values[i]

                                    # Fill empty pixels with defined Gaussian noise (mean=0, std=1)
                                    mask = Image == 0

                                    # Need to add noise for proper prediction
                                    noise_pix = np.random.normal(background, noise, Image.shape) 
                                    Image[mask] = noise_pix[mask]

                                    # Normalize each sample by its own mean and std
                                    Images_norm = IndividualNormalize(Image)

                                    upsampling_factor = 10
                                    # Upsampling using a simple nearest neighbor interp.
                                    Image_upsampled = np.kron(Images_norm, np.ones((1, upsampling_factor, upsampling_factor)))
                                    
                                    # Convert the numpy array to a PyTorch tensor and move it to the CUDA device
                                    Images_norm_tensor = torch.as_tensor(Image_upsampled).float().unsqueeze(0).to(device)
                                    
                                    predicted_density = model(Images_norm_tensor)

                                    predicted_density_cpu = predicted_density.cpu().detach().numpy()
                                    predicted_density_cpu = predicted_density_cpu.squeeze()

                                    # Threshold negative values
                                    predicted_density_cpu[predicted_density_cpu < 0] = 0

                                    SRmask, SRArea = simple_mask_cutoff(pred_img=predicted_density_cpu, min_threshold=0.1)

                                    guess_D = estimate_D(SRArea,coeff_a=1162.0)

                                    if guess_D > 1:                           
                                        # Use weighted centroid
                                        # Estimate localization
                                        # Calculate center of mass using intensity values as weights
                                        y0, x0 = regionprops(SRmask.astype(int), intensity_image=predicted_density_cpu)[0].centroid_weighted
                                        x0 = x0/10-0.5 # revert back to 110 nm pixel size pixel coordinate
                                        y0 = y0/10-0.5 

                                        # Convert to the same pixel coordinate as x_values and y_values from raw image coordinate
                                        x1 = x0 - x_offset # matlab starting from 1 pixel coordinate do not need to -1
                                        y1 = y0 - y_offset 

                                        # Add a new row to the DataFrame
                                        df_temp = pd.DataFrame({'Frame': iFrame, 'PSF_idx': iPSF, 'Xpos': x1, 'Ypos': y1,
                                                        'TotalPhoton': tab_psf_fitresult.TotalPhoton[irow], 'Intensity': tab_psf_fitresult.Intensity[irow],
                                                        'Background': tab_psf_fitresult.Background[irow],'EllipticityIdx': tab_psf_fitresult.EllipticityIdx[irow],'Angle': tab_psf_fitresult.Angle[irow], 
                                                        'SNR': tab_psf_fitresult.SNR[irow],'COV': tab_psf_fitresult.COV[irow],'UNetArea': UNetArea, 'SRArea': SRArea, 'D': guess_D}, index=[0])
                                        df_temp.to_csv(path.join(csv_save_path, sr_fileName), header=False, index=False, mode='a')
                                    else:
                                        # Add a new row to the DataFrame
                                        df_temp = pd.DataFrame({'Frame': iFrame, 'PSF_idx': iPSF, 'Xpos': tab_psf_fitresult.Xpos[irow], 'Ypos': tab_psf_fitresult.Ypos[irow],
                                                        'TotalPhoton': tab_psf_fitresult.TotalPhoton[irow], 'Intensity': tab_psf_fitresult.Intensity[irow],
                                                        'Background': tab_psf_fitresult.Background[irow],'EllipticityIdx': tab_psf_fitresult.EllipticityIdx[irow],'Angle': tab_psf_fitresult.Angle[irow], 
                                                        'SNR': tab_psf_fitresult.SNR[irow],'COV': tab_psf_fitresult.COV[irow],'UNetArea': UNetArea, 'SRArea': SRArea, 'D': guess_D}, index=[0])
                                        df_temp.to_csv(path.join(csv_save_path, sr_fileName), header=False, index=False, mode='a')                             
                                irow += 1
                        iFrame += 1



                
                # Handle the remaining frames (if any)  
                remaining_frames = ND2_img.shape[0] % slice_interval  
                if remaining_frames > 0:  
                    start_index = num_full_slices * slice_interval  
                    frame_slice = ND2_img[start_index:, :, :]  
                    # Perform your analysis on the remaining frames  
                    print(f"Shape of remaining frame slice: {frame_slice.shape}")  
                    bk_img = estimate_background_image(frame_slice,device)
                    torch.cuda.empty_cache()
                    # bk_img = temp_bk_img                

                    for isubFrame in tqdm(range(remaining_frames)):
            
                        xyI_dataset = cell_PSF[cell_PSF['xyI'][0][isubFrame+start_index]]  # last index is frame number
                            
                        total_PSF_num = xyI_dataset[0].size
                        if xyI_dataset.len() == 1: # this frame contains PSF, if 2 means empty frame       

                            for iPSF in range(total_PSF_num):
                            
                                referenced_data = matfile[xyI_dataset[0, iPSF]]  # last index is the PSF index in the current frame
                                xyI_array = np.array(referenced_data)  # Convert to numpy array
                                UNetArea = xyI_array.shape[1] # Store the U-Net area

                                # Extract x, y, and intensity values
                                x_values = xyI_array[0]
                                y_values = xyI_array[1]
                                intensity_values = xyI_array[2]
                                
                                center_column = (x_values.min()+x_values.max())/2
                                center_row = (y_values.min()+y_values.max())/2
                                    
                                # Place the original array into the center of the new array
                                x_offset = int(32 / 2) - round(center_column)
                                y_offset = int(32 / 2) - round(center_row)

                                # Create a 32x32 array
                                Image = np.zeros((32, 32))

                                # Check if out-of-Image dimension, discard if it is too big, most like UNet mask cover more than 2 molecules
                                if (y_values+y_offset).min() >= 0 and (y_values+y_offset).max() < 32 and (x_values+x_offset).min() >= 0 and (x_values+x_offset).max() < 32:

                                    # Get background and noise from bk_img
                                    background = np.mean(bk_img[isubFrame,y_values.astype(int)-1,x_values.astype(int)-1])
                                    noise = np.std(bk_img[isubFrame,y_values.astype(int)-1,x_values.astype(int)-1])
                                    
                                    for i in range(len(x_values)):
                                        Image[int(y_values[i]) + y_offset,int(x_values[i]) + x_offset] = intensity_values[i]

                                    # Fill empty pixels with defined Gaussian noise (mean=0, std=1)
                                    mask = Image == 0

                                    # Need to add noise for proper prediction
                                    noise_pix = np.random.normal(background, noise, Image.shape) 
                                    Image[mask] = noise_pix[mask]

                                    # Normalize each sample by its own mean and std
                                    Images_norm = IndividualNormalize(Image)

                                    upsampling_factor = 10
                                    # Upsampling using a simple nearest neighbor interp.
                                    Image_upsampled = np.kron(Images_norm, np.ones((1, upsampling_factor, upsampling_factor)))
                                    
                                    # Convert the numpy array to a PyTorch tensor and move it to the CUDA device
                                    Images_norm_tensor = torch.as_tensor(Image_upsampled).float().unsqueeze(0).to(device)
                                    
                                    predicted_density = model(Images_norm_tensor)

                                    predicted_density_cpu = predicted_density.cpu().detach().numpy()
                                    predicted_density_cpu = predicted_density_cpu.squeeze()

                                    # Threshold negative values
                                    predicted_density_cpu[predicted_density_cpu < 0] = 0

                                    SRmask, SRArea = simple_mask_cutoff(pred_img=predicted_density_cpu, min_threshold=0.1)

                                    guess_D = estimate_D(SRArea,coeff_a=1162.0)

                                    if guess_D > 1:                           
                                        # Use weighted centroid
                                        # Estimate localization
                                        # Calculate center of mass using intensity values as weights
                                        y0, x0 = regionprops(SRmask.astype(int), intensity_image=predicted_density_cpu)[0].centroid_weighted
                                        x0 = x0/10-0.5 # revert back to 110 nm pixel size pixel coordinate
                                        y0 = y0/10-0.5 

                                        # Convert to the same pixel coordinate as x_values and y_values from raw image coordinate
                                        x1 = x0 - x_offset # matlab starting from 1 pixel coordinate do not need to -1
                                        y1 = y0 - y_offset 

                                        # Add a new row to the DataFrame
                                        df_temp = pd.DataFrame({'Frame': iFrame, 'PSF_idx': iPSF, 'Xpos': x1, 'Ypos': y1,
                                                        'TotalPhoton': tab_psf_fitresult.TotalPhoton[irow], 'Intensity': tab_psf_fitresult.Intensity[irow],
                                                        'Background': tab_psf_fitresult.Background[irow],'EllipticityIdx': tab_psf_fitresult.EllipticityIdx[irow],'Angle': tab_psf_fitresult.Angle[irow], 
                                                        'SNR': tab_psf_fitresult.SNR[irow],'COV': tab_psf_fitresult.COV[irow],'UNetArea': UNetArea, 'SRArea': SRArea, 'D': guess_D}, index=[0])
                                        df_temp.to_csv(path.join(csv_save_path, sr_fileName), header=False, index=False, mode='a')
                                    else:
                                        # Add a new row to the DataFrame
                                        df_temp = pd.DataFrame({'Frame': iFrame, 'PSF_idx': iPSF, 'Xpos': tab_psf_fitresult.Xpos[irow], 'Ypos': tab_psf_fitresult.Ypos[irow],
                                                        'TotalPhoton': tab_psf_fitresult.TotalPhoton[irow], 'Intensity': tab_psf_fitresult.Intensity[irow],
                                                        'Background': tab_psf_fitresult.Background[irow],'EllipticityIdx': tab_psf_fitresult.EllipticityIdx[irow],'Angle': tab_psf_fitresult.Angle[irow], 
                                                        'SNR': tab_psf_fitresult.SNR[irow],'COV': tab_psf_fitresult.COV[irow],'UNetArea': UNetArea, 'SRArea': SRArea, 'D': guess_D}, index=[0])
                                        df_temp.to_csv(path.join(csv_save_path, sr_fileName), header=False, index=False, mode='a')                                   
                                irow += 1
                        iFrame += 1    
                ND2_index_adjusted += 1
              

    # # write results to a csv
    # tif_save_path = path.join(rootDir,nowDir)
    # df.to_csv(path.join(tif_save_path, sr_fileName), index=False)