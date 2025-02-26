# SMLMM

Single-Molecule Localization and Mobility Microscopy

# Installation

SMLMM has the following depencies:

Python 3.7

Fiji distribution of ImageJ

MATLAB 2021 or higher version

Hardware: NVIDIA GPU supporting CUDA 11.7, we use NVIDIA Quadro P4000.

We recommend create a new conda environment with Python version 3.7, using the command below. After activating the conda environment, cd to the directory created from the Github clone command, and use below command to install the dependent packages from requirements.txt.

```
pip install -r requirements.txt
```

For deep learning training and inference, we use CUDA 11.7 with PyTorch 1.13.0. You can install PyTorch using the following command:

```
# CUDA 11.7
conda install pytorch==1.13.0 torchvision==0.14.0 torchaudio==0.13.0 pytorch-cuda=11.7 -c pytorch -c nvidia
```

# Train your own models

You can download our pretrained models: the **UNet model** ([checkpoint_epoch13.pth](https://1drv.ms/u/s!ApPGm6eczDkKgUnLgOphRpAwutgG?e=VMcofl)) and the **Deep-SnapTrack model** ([MBX_20231220_110nmPix_rep2_epoch9.pth](https://1drv.ms/u/s!ApPGm6eczDkKgUqF7B__WWXfqQRU?e=IYTqNz)) from the provided OneDrive links. Please note that these models are trained under specific camera setups and experimental conditions. Depending on your own setup and data, you may need to retrain these models to achieve optimal performance for your application. Below is a guide to help you train your own models if needed.

## Step 1: simSnapshot

A pipeline to convert ground truth trajectories into simulated images for training process. The tracks are generated using a pubished method called simSPT ([Tjian - Darzacq lab / simSPT · GitLab](https://gitlab.com/tjian-darzacq-lab/simSPT)). We generated a set trajectories with logarithmically increasing diffusion coefficient through [simSPT_script.txt](https://1drv.ms/t/s!ApPGm6eczDkKgUyfGX_uYKiVVq7R?e=s2iwdF) in the command line. For convenience, we simulated long tracks and segmented them into shorter segments suitable for the target exposure time. Notably, we modified the `betaUnit` parameter in the original `simSPT.c` file from 0 to 1 to count track lifetime in seconds.

Next, we generated training data pairs for the UNet model. First, we converted the x-y coordinates of the tracks into pixel positions and overlaid four tracks onto a single 32×32 image. The images were then binarized and convolved with a weight map to generate ground truth masks. To mimic real camera acquisition, we added Gaussian white noise and Poisson shot noise to the images and paired them with ground truth masks. This process can be executed using the provided script:

`step0_sim_4molecules_forUNet.m`

For generating training data pairs for DeepSnapTrack, we converted the x-y coordinates of the tracks into pixel positions with 1/10 of the size used for UNet. Each track was overlaid onto a separate 320×320 image. The ground truth masks were defined as the top 95% of the total signal volume. Simulated images were then generated using the same noise augmentation method as in UNet, but resized to 32×32 pixels. This process can be executed using the provided script:

`step0_sim_1molecule_forDeepSnapTrack.m`

## Step 2: UNet training

After preparing single molecule snapshot images and corresponding masks according to **Step 0: Sim Snapshot**. Specify your path to data pairs in `step1-0_train_multiGPU.sh` by defining:

```
--dir_img  /path/to/snapshot/images \
--dir_mask  /path/to/masks \
```

And define a data path where you save your checkpoint files:

```
--dir_checkpoint /path/to/save/checkpoints \
```

Run the `.sh` file to generate customized UNet models.

## Step 3: Deep-SnapTrack training

In order to train Deep-SnapTrack with your own data, you need to prepare simulated single molecule snapshot images and corresponding ground truth trajectory maps. Specify your path to data pairs in `step3-0_ptorch_train.sh` by defining:

```
--dir_img  /path/to/snapshot/images \
--dir_mask  /path/to/trajectory/maps \
```

And define a data path where you save your checkpoint files:

```
--dir_checkpoint /path/to/save/checkpoints \
```

Run the `.sh` file to generate customized Deep-SnapTrack models.

# SMLMM data analysis

## Step 1: UNet Segmentation

To do segmentation with UNet, specify your data path in the corresponding section of `step1-1_ND2batch_prediction.sh`:

    --model /path/to/your/model/checkpoint_epoch13.pth \
    --input /path/to/your/data \
    --output /path/to/save/masks

 The output will be a `.tif` file corresponding to your raw data named `UNet_mask_MBX_20240620_epoch13_Ch1.tif`

## Step 2: Single Molecule Selection

### Do ThunderSTORM Localizations

We use ThunderSTORM (https://github.com/zitmen/thunderstorm) to assist in selecting qualified snapshot from UNet result. To do this, you should first make sure you have download ThunderSTORM plugin in your imagej. Then, choose `Run Macro` and select `step2-1_IJmacro_ThunderSTORM.ijm`, which will open an interface for you to load raw images and define output path for localization result. To note, the previous step will create a sub folder named ThunderSTORM under "mask" folder:

```
/path/to/save/masks/ThunderSTORM
```

Your output path are suggested to be defined as this.

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

Next you have to specify your data path:

```matlab
% >>>>>>>>>>>>>>>>>>>> MOTION BLUR DETECTION PARAMETERS >>>>>>>>>>>>>>>>>>>> %
% directory of your raw image data (ND2 file)
input_path = '/path/to/your/data/';

% the parent folder of motion blur detection and analysis
output_path = '/path/to/save/masks/';

% defined name of your result for each cell
input_rawND2_prefix = {...
    'name_of_result',
     };
% <<<<<<<<<<<<<<<<<<<< MOTION BLUR DETECTION PARAMETERS <<<<<<<<<<<<<<<<<<<<< %
```

After finishing these settings, run `step2-2_script_snapshot_detection.m`. The output will be a file named `Blurdata_UNet_mask_MBX_20240620_epoch13_Ch1.mat`. Note that all names have to be defined in `input_rawND2_prefix` if your raw data contains multiple clusters.

## Step 3: Deep-SnapTrack

To do prediction with Deep-SnapTrack, specify the data path in the corresponding section of `step3-1_batch_prediction.py`

```python
# >>>>>>>>>>>>>>>>>>>>>>>>>> input >>>>>>>>>>>>>>>>>>>>>>>>>> #
# directory of your DeepSnapTrack model
weights_file = '/path/to/your/model/MBX_20231220_110nmPix_rep2_epoch9.pth'

# parent directory where you save UNet motion blur extracted mat folder
rootDir = '/path/to/save/masks'

# directory of UNet motion blur extracted mat folder
dataDir = [   
    'name_of_folder_with_.mat_file',
    ]

UNet_model = 'UNet_mask_MBX_20240620_2035_epoch20_Ch1'
blurmat_file_prefix = 'Blurdata_'+UNet_model
fitresult_file = 'Fitresult_'+UNet_model+'.csv'
sr_fileName = UNet_model+'_SR_pred_v3.csv'

# raw image files, keep order same as dataDir
ND2File = [    
    '/path/to/your/data',
    ]

# <<<<<<<<<<<<<<<<<<<<<<<<<<<< input <<<<<<<<<<<<<<<<<<<<<<<<<<<< #
```

Then run the following command to generate final results

```
python step3-1_batch_prediction.py
```

The output will be a CSV file named `UNet_mask_MBX_20240620_epoch13_Ch1_SR_pred_v3.csv`, saved in your input `dataDir`, along with your previously generated MAT file named `Blurdata_UNet_mask_MBX_20240620_epoch13_Ch1.mat`.

### Step 4: MPALM rendering

 You can get the final visualization result by uploading the two files generated from the last step` UNet_mask_MBX_20240620_epoch13_Ch1_SR_pred_v3.csv`  and `Blurdata_UNet_mask_MBX_20240620_epoch13_Ch1.mat `to the MATLAB app provided under ./step4_MPALM_rendering/step4_main_mobilityPALM.mlapp. For the user guide of this app, please see ./step4_MPALM_rendering/Users Guide.docx.
