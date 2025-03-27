# SMLDM 

Single Molecule Localization and Diffusivity Microscopy

# Installation

SMLDM has the following depencies:

Python 3.7

Fiji distribution of ImageJ

MATLAB 2021 or higher version

Hardware: NVIDIA GPU supporting CUDA 11.7, we use NVIDIA Quadro P4000.

### Seting up a conda environment

We recommend create a new conda environment with Python version 3.7. After activating the conda environment, cd to the directory created from the Github clone command, and use below command to install the dependent packages from requirements.txt.

```bash
pip install -r requirements.txt
```

For deep learning training and inference, we use CUDA 11.7 with PyTorch 1.13.0. In the activated conda environment, you can install PyTorch using the following command:

```bash
# CUDA 11.7
conda install pytorch==1.13.0 torchvision==0.14.0 torchaudio==0.13.0 pytorch-cuda=11.7 -c pytorch -c nvidia
```

### Setting Up CUDA Environment Path

Please download the CUDA Toolkit (11.7) from Nvidia website https://developer.nvidia.com/cuda-11-7-0-download-archive and install it in your system according to the website instruction.

To enable PyTorch to run on the GPU, you need to add the CUDA toolkit to the system environment path. In linux system, follow these steps to create activate and deactivate scripts for managing the CUDA environment variables:

1. **Create the Activation Script**
   - Create the directory for activation scripts:
     ```bash
     mkdir -p $CONDA_PREFIX/etc/conda/activate.d
     ```
   - Create a file named `env_vars.sh` in the activation directory with the following content:
     ```bash
     #!/bin/sh
     # Environment variable for CUDA 11.7
     export PATH=/usr/local/cuda-11.7/bin${PATH:+:${PATH}}
     export LD_LIBRARY_PATH=/usr/local/cuda-11.7/lib64${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}
     ```

2. **Create the Deactivation Script**
   - Create the directory for deactivation scripts:
     ```bash
     mkdir -p $CONDA_PREFIX/etc/conda/deactivate.d
     ```
   - Create a file named `env_vars.sh` in the deactivation directory with the following content:
     ```bash
     #!/bin/sh
     # Unset environment variable for CUDA 11.7
     export PATH=${PATH#/usr/local/cuda-11.7/bin:}
     unset LD_LIBRARY_PATH
     ```

These scripts will automatically set and unset the CUDA environment variables when you activate and deactivate your conda environment, respectively.


# Train your own models

You can download our provided the pretrained U-Net model (`checkpoint_UNet_epoch20.pth`) and the pretrained Deep-SnapTrack model (`MBX_20231220_110nmPix_rep2_epoch9.pth`) from the zenodo repository https://zenodo.org/records/15089354. Please note that these models are trained under specific camera setups (camera pixel size 110 nm) and experimental conditions (exposure time 30 ms). Depending on your own setup and data, you recommend measure your microscopic point-spread-function using the fluorescent beads and generate simulated training dataset as explained below, and to retrain these models to achieve optimal performance for your application. 

## Step 1: simSnapshot

We developed a pipeline to convert the simulated Brownian diffusing trajectories into simulated molecule snapshot images together with its mask and track image as datasets for downstream model training. The molecule trajectories were generated using a pubished method called simSPT ([Tjian - Darzacq lab / simSPT · GitLab](https://gitlab.com/tjian-darzacq-lab/simSPT)). We generated a set of trajectories with logarithmically increasing diffusion coefficients using the command recorded in `simSPT_script.txt`. For convenience, we simulated long tracks with 1-ms time interval and segmented them into shorter segments matching the desired exposure time. Notably, we modified the `betaUnit` parameter in the original `simSPT.c` file from 0 to 1 to count track lifetime in seconds.

We provided our simulated molecule trajectories under the zenodo repository https://zenodo.org/records/15089354.

Next, we generated training datasets for the U-Net segmentation model. First, we placed four molecule trajectories into a single 32×32 image with a pixel size of 110 nm. The positions of the molecule trajectories were convolved with a Gaussian model of the microscopic point-spread function (PSF). The resulting matrix was then binned, rescaled to the desired signal-to-noise ratio (SNR), and augmented with both Gaussian white noise and Poisson shot noise to mimic experimental molecule snapshots. The mask for each molecule was defined as the region encompassing the top 95% of the total intensity signal. The simulated snapshot and corresponding mask images were paired to serve as the ground truth for U-Net training. This process can be executed using the following script:

`step0_sim_4molecules_forUNet.m`

To generate training data pairs for Deep-SnapTrack, we placed a single molecule trajectory into a 32×32 image with a pixel size of 110 nm and generated its snapshot image using the method described above. The molecule trajectories were then convolved with a 7×7 Gaussian kernel (σ = 1 pixel-width), and this convolution was evaluated at an 11-nm pixel width to produce a higher-resolution 320×320 track image. The final training dataset for Deep-SnapTrack consisted of paired simulated snapshot images and their corresponding track images. This process can be executed using the following script:

`step0_sim_1molecule_forDeepSnapTrack.m`

## Step 2: U-Net training

After preparing single molecule snapshot images and corresponding masks according to **Step 0: Sim Snapshot**. Specify your path to data pairs in `step1-0_train_multiGPU.sh` by defining:

```
--dir_img  /path/to/save/simulated_dataset_forUNet/images/ \
--dir_mask  /path/to/save/simulated_dataset_forUNet/masks/ \
```

And define a data path where you save your checkpoint files:

```
--dir_checkpoint ./checkpoints/ \
```

Run the `.sh` file to generate customized U-Net models.

## Step 3: Deep-SnapTrack training

In order to train Deep-SnapTrack with your own data, you need to prepare simulated single molecule snapshot images and corresponding ground truth trajectory maps. Specify your path to data pairs in `step3-0_ptorch_train.sh` by defining:

```
--dir_img  /path/to/save/simulated_dataset_forDeepSnapTrack/imgs \
--dir_mask  /path/to/save/simulated_dataset_forDeepSnapTrack/trackHeatmap \
```

And define a data path where you save your checkpoint files:

```
--dir_checkpoint ./checkpoints/ \
```

Run the `.sh` file to generate customized Deep-SnapTrack models.

# SMLDM data analysis

We will use one ND2 image captured for Paxillin-Halo as the example to demonstrate the analysis pipeline. You can download this image file from the the zenodo repository https://zenodo.org/records/15089354.

## Step 1: U-Net Segmentation

To do segmentation with U-Net, specify your data path in the corresponding section of `step1-1_ND2batch_prediction.sh`:

    --model /path/to/your/model/checkpoint_UNet_epoch20.pth \
    --input /path/to/your/data/20240712_Clust01_U2OS_Paxillin_30p5ms_2kframe_001.nd2 \
    --output /path/to/save/results/

 The output will be a `.tif` file corresponding to your raw data named `UNet_mask_MBX_20240620_epoch20_Ch1.tif`

## Step 2: Single Molecule Selection

### Do ThunderSTORM Localizations

We use ThunderSTORM (https://github.com/zitmen/thunderstorm) to assist in selecting qualified snapshot from U-Net result, i.e. only one ThunderSTORM localization in a U-Net mask is kept as single-molecule snapshot for downstream analysis. To do this, please first install [ThunderSTORM](zitmen.github.io/thunderstorm/) plugin in your ImageJ. Then, choose `Run Macro` and select `step2-1_IJmacro_ThunderSTORM.ijm`, which will open an interface for you to load raw images and define output path for localization result. Please create a new folder called ThunderSTORM under the parent folder of previous output folder, and save the result of ThunderSTORM in this folder:

```
/path/to/save/results/ThunderSTORM
```
After finish the ThunderSTORM step, move to `step2-2_script_snapshot_detection.m` to use ThunderSTORM localization filter the qualified single-molecule snapshots.

### Select Qualified Snapshots

If your data was generated with bulk images captured after each illumination sequence, for example after every 2000 frames, you can set `has_bulk_image` true and manually select an ROI. Otherwise you can set `draw_fullROI` true and no masks will be eliminated in this step.

```matlab
% >>>>>>>>>>>>>>>>>>>> NUCLEUS SELECTION >>>>>>>>>>>>>>> %
has_bulk_imgs = false;
skip_drawROI = false; % set true if already done
draw_fullROI = true;
WideField_subfolder = 'BF';
% <<<<<<<<<<<<<<<<<<<< NUCLEUS SELECTION <<<<<<<<<<<<<<< %
```

Next, specify your data path:

```matlab
% >>>>>>>>>>>>>>>>>>>> MOTION BLUR DETECTION PARAMETERS >>>>>>>>>>>>>>>>>>>> %
% directory of your raw image data (ND2 file)
input_path = '/path/to/your/data/';

% the parent folder of motion blur detection and analysis
output_path = '/path/to/save/results/';

% defined name of your result for each cell
input_rawND2_prefix = {...
    '20240712_Clust01_U2OS_Paxillin',
     };
% <<<<<<<<<<<<<<<<<<<< MOTION BLUR DETECTION PARAMETERS <<<<<<<<<<<<<<<<<<<<< %
```

After finishing these settings, run `step2-2_script_snapshot_detection.m`. The output will be a file named `Blurdata_UNet_mask_MBX_20240620_epoch20_Ch1.mat`. Note that all names have to be defined in `input_rawND2_prefix` if your raw data contains multiple image sequence files.

## Step 3: Deep-SnapTrack

To do prediction with Deep-SnapTrack, specify the data path in the corresponding section of `step3-1_batch_prediction.py`

```python
# >>>>>>>>>>>>>>>>>>>>>>>>>> input >>>>>>>>>>>>>>>>>>>>>>>>>> #
# directory of your DeepSnapTrack model
weights_file = './checkpoints/MBX_20231220_110nmPix_rep2_epoch9.pth'

# parent directory where you save U-Net motion blur extracted mat folder
rootDir = '/path/to/save/results'

# directory of U-Net motion blur extracted mat folder
dataDir = [   
    '20240712_Clust01_U2OS_Paxillin_Cell01',
    ]

UNet_model = 'UNet_mask_MBX_20240620_epoch20_Ch1'
blurmat_file_prefix = 'Blurdata_'+UNet_model
fitresult_file = 'Fitresult_'+UNet_model+'.csv'
sr_fileName = UNet_model+'_SR_pred_v3.csv'

# raw image files, keep order same as dataDir
ND2File = [    
    '/path/to/your/data/20240712_Clust01_U2OS_Paxillin_30p5ms_2kframe_001.nd2',
    ]

# <<<<<<<<<<<<<<<<<<<<<<<<<<<< input <<<<<<<<<<<<<<<<<<<<<<<<<<<< #
```

Then run the following command to generate final results

```
python step3-1_batch_prediction.py
```

The output will be a CSV file named `UNet_mask_MBX_20240620_epoch20_Ch1_SR_pred_v3.csv`, saved in your input `dataDir`, along with your previously generated MAT file named `Blurdata_UNet_mask_MBX_20240620_epoch20_Ch1.mat`.

The csv file contains:

| Column            | Description                                                                 | Unit                     |
|-----------------|-----------------------------------------------------------------------------|--------------------------|
| Frame           | Frame number                                                                | -                       |
| PSF_idx         | Index of snapshot                                                           | -                       |
| Xpos            | X coordinate of molecule                                 | Image pixel              |
| Ypos            | Y coordinate of molecule                                 | Image pixel              |
| TotalPhoton     | Estimated total photon number (unit: photons)                               | Photons                  |
| Intensity       | Intensity value from elliptical Gaussian fitting                | a.u.                     |
| Background      | Background value from elliptical Gaussian fitting                | a.u.                     |
| EllipticityIdx  | Ellipticity index (no longer used, please ignore)                            | -                       |
| Angle           | Angle of ellipticity (no longer used, please ignore)                        | -                       |
| SNR             | Signal-to-noise ratio                                       | Decibel                  |
| COV             | Covariance (another measurement for SNR)                                    | -                       |
| UNetArea        | Molecule snapshot area from U-Net masks                   | Image pixel              |
| SRArea          | Molecule pseudo-track area, can convert to diffusion coeffcient                     | 1/10 of image pixel      |
| D          | Molecule diffusion coefficient                      | um^2/s      |

The mat file contains:
| Variable Name       | Description                                                                 | Variable Type |
|---------------------|-----------------------------------------------------------------------------|---------------|
| UNet_model_folder   | Path of U-Net masks                                                         | cell          |
| cell_PSF            | Stores the pixel coordinate and pixel intensity of detected snapshots       | struct        |
| filter              | Filter parameters defined in step2-2_script_snapshot_detection.m          | struct        |
| impars              | Image acquisition parameters defined in step2-2_script_snapshot_detection.m | struct        |
| input_path          | Raw image path                                                              | char          |
| input_rawND2_prefix | The prefix of image used for naming files                                   | char          |
| output_path         | Path to save result                                                         | char          |
### Step 4: MPALM rendering

 You can get the final visualization result by uploading the two files generated from the last step` UNet_mask_MBX_20240620_epoch20_Ch1_SR_pred_v3.csv`  and `Blurdata_UNet_mask_MBX_20240620_epoch20_Ch1.mat `to the MATLAB app provided under ./step4_MPALM_rendering/step4_main_mobilityPALM.mlapp. For the user guide of this app, please see ./step4_MPALM_rendering/Users Guide.docx.
