% This is a script to prepare input for pooled MPALM analysis
%% PUT YOU INPUT BELOW
%% ---------------- 1xNLS/6xNLS/H2B ------------------%
clearvars data_struct FinalResults
% general parameters
FitParams.DoSingleFit = false;
FitParams.DoSinglePlot = false;
FitParams.DoMergeFit = false;
FitParams.DoMergePlot = true;
% hist_edges = linspace(0,4,100);

% load sample
iSamp = 1;
rootPath = 'E:\Submission\2023_MotionBlur\PA646MBX_H2B_1x6xNLS_beads\1xNLS\UNet_mask_MBX_20231203_1016_epoch5';
name_pattern = "30p5ms_20kframe";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Result_UNet_mask_MBX_20231203_1016_epoch5.mat";
data_struct(iSamp).SampleName = "1xNLS";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20231203_1016_epoch5_SR_pred_v2.csv";%"SR_pred.csv";

iSamp = 2;
rootPath = 'E:\Submission\2023_MotionBlur\PA646MBX_H2B_1x6xNLS_beads\6xNLS\UNet_mask_MBX_20231203_1016_epoch5';
name_pattern = "30p5ms_5kframe";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Result_UNet_mask_MBX_20231203_1016_epoch5.mat";
data_struct(iSamp).SampleName = "6xNLS";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20231203_1016_epoch5_SR_pred_v2.csv";%"SR_pred.csv";

iSamp = 3;
rootPath = 'E:\Submission\2023_MotionBlur\PA646MBX_H2B_1x6xNLS_beads\H2B\UNet_mask_MBX_20231203_1016_epoch5';
name_pattern = "30p5ms_5kframe";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Result_UNet_mask_MBX_20231203_1016_epoch5.mat";
data_struct(iSamp).SampleName = "H2B";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20231203_1016_epoch5_SR_pred_v2.csv";%"SR_pred.csv";

iSamp = 4;
rootPath = 'E:\Submission\2023_MotionBlur\PA646MBX_H2B_1x6xNLS_beads\beads';
name_pattern = "30p5ms_5kframe_01";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Result_UNet_mask_MBX_20230712_1000_epoch12.mat";
data_struct(iSamp).SampleName = "beads";
data_struct(iSamp).SRFile = "SRpred_singleBeads_MBX_20231220_110nmPix_rep2_epoch9.csv";

% iSamp = 5;
% rootPath = 'E:\PROCESS-SPT\2024\20240222_U2OS_H2B_1xNLS_PA646_MPALM';
% name_pattern = "1xNLS_Cell";
% temp_file = dir(fullfile(rootPath));
% temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern);
% data_struct(iSamp).workspaces = temp_file(idx);
% data_struct(iSamp).UNetMat = "Result_Ch1_UNet_mask_MBX_20240103_1400_epoch19_Ch1.mat";
% data_struct(iSamp).SampleName = "20240222 1xNLS";
% 
% iSamp = 6;
% rootPath = 'E:\PROCESS-SPT\2024\20240222_U2OS_H2B_1xNLS_PA646_MPALM';
% name_pattern1 = "Cell01"; name_pattern2 = "H2B";name_pattern3 = "Clust01";
% temp_file = dir(fullfile(rootPath));
% temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2) & contains(temp_Filenames,name_pattern3);
% data_struct(iSamp).workspaces = temp_file(idx);
% data_struct(iSamp).UNetMat = "Result_Ch1_UNet_mask_MBX_20240103_1400_epoch19_Ch1.mat";
% data_struct(iSamp).SampleName = "20240222 H2B";
% 
% iSamp = 7;
% rootPath = 'E:\PROCESS-SPT\2024\20240125_U2OS_HaloRPB1_HMSiR_PA646\HMSiR10nM_TMR';
% name_pattern = "Cell";
% temp_file = dir(fullfile(rootPath));
% temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern);
% data_struct(iSamp).workspaces = temp_file(idx);
% data_struct(iSamp).UNetMat = "Result_Ch1_UNet_mask_MBX_20240103_1400_epoch19_Ch1.mat";
% data_struct(iSamp).SampleName = "20240125 RPB1";

%%
clearvars data_struct FinalResults
% general parameters
FitParams.DoSingleFit = false;
FitParams.DoSinglePlot = true;
FitParams.DoMergeFit = true;
FitParams.DoMergePlot = true;
% hist_edges = linspace(0,4,100);

iSamp = 1;
rootPath = 'E:\PROCESS-SPT\2024\20240320_U2OS_EF1a_H2B-Halo_PA646_DRB_FLA\30p5ms';
name_pattern1 = "DMSO"; name_pattern2 = "Cell";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20231203_1016_epoch5_Ch1.mat";
data_struct(iSamp).SampleName = "DMSO";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20231203_1016_epoch5_Ch1_SR_pred.csv";

iSamp = 2;
rootPath = 'E:\PROCESS-SPT\2023\20231013_U2OS_H2BOnly_HP1-H2B_mPALM\EF1a-H2B';
name_pattern1 = "H2B"; name_pattern2 = "Cell";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20231203_1016_epoch5.mat";
data_struct(iSamp).SampleName = "2023 H2B";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20231203_1016_epoch5_SR_pred.csv";

iSamp = 2;
rootPath = 'E:\PROCESS-SPT\2024\20240320_U2OS_EF1a_H2B-Halo_PA646_DRB_FLA\30p5ms';
name_pattern1 = "10uMFLA"; name_pattern2 = "Cell";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20231203_1016_epoch5_Ch1.mat";
data_struct(iSamp).SampleName = "10uM FLA";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20231203_1016_epoch5_Ch1_SR_pred.csv";

iSamp = 3;
rootPath = 'E:\PROCESS-SPT\2024\20240320_U2OS_EF1a_H2B-Halo_PA646_DRB_FLA\30p5ms';
name_pattern1 = "100uMDRB"; name_pattern2 = "Cell";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20231203_1016_epoch5_Ch1.mat";
data_struct(iSamp).SampleName = "100uM DRB";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20231203_1016_epoch5_Ch1_SR_pred.csv";

%% ---------------- RPB1 DMOS/THZ1 ------------------%

% general parameters
FitParams.DoSingleFit = false;
FitParams.DoSinglePlot = false;
FitParams.DoMergeFit = False;
FitParams.DoMergePlot = true;
% hist_edges = 0:0.05:2;
% hist_edges = 20:5:200;

% load sample
clearvars data_struct FinalResults

iSamp = 1;
rootPath = 'E:\PROCESS-SPT\2023\20231115_U2OS_RPB1_fastSPT_mPALM_FDAligned';
name_patternA = "Cell0"; name_patternB = "DMSO";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_patternA) & contains(temp_Filenames,name_patternB);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Result_v2_UNet_mask_MBX_20230903_2330_epoch5.mat";
data_struct(iSamp).SampleName = "DMSO";

iSamp = 2;
rootPath = 'E:\PROCESS-SPT\2023\20231115_U2OS_RPB1_fastSPT_mPALM_FDAligned';
name_patternA = "Cell0"; name_patternB = "1uMTHZ1";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_patternA) & contains(temp_Filenames,name_patternB);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Result_v2_UNet_mask_MBX_20230903_2330_epoch5.mat";
data_struct(iSamp).SampleName = "1uMTHZ1";

iSamp = 1;
rootPath = 'E:\PROCESS-SPT\2024\20240118_U2OS_HaloRPB1_JF549_JFX640_MPALM\HTL_JFX650_TMR';
name_patternA = "Cell0"; 
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_patternA);
data_struct(iSamp).workspaces = temp_file(idx); data_struct(iSamp).workspaces = data_struct(iSamp).workspaces(5:20);
data_struct(iSamp).UNetMat = "Result_Ch2_UNet_mask_MBX_20231212_1300_epoch3_Ch2.mat";
data_struct(iSamp).SampleName = "JFX650_25nM";

%% High and low localization density for H2B and 1xNLS
clearvars data_struct FinalResults
% general parameters
FitParams.DoSingleFit = true;
FitParams.DoSinglePlot = false;
FitParams.DoMergeFit = true;
FitParams.DoMergePlot = true;

% load sample
iSamp = 1;
rootPath = 'F:\Submission\2024_MotionBlur\PA646MBX_H2B_1x6xNLS_beads\H2B\UNet_mask_MBX_20231203_1016_epoch5';
name_pattern = "30p5ms_5kframe";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Result_UNet_mask_MBX_20231203_1016_epoch5.mat";
data_struct(iSamp).SampleName = "sparse H2B";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20231203_1016_epoch5_SR_pred_v2.csv";%"SR_pred.csv";
FinalResults{iSamp}.locPerFrame = [];

iSamp = 2;
rootPath = 'F:\PROCESS-SPT\2024\20240320_U2OS_EF1a_H2B-Halo_PA646_DRB_FLA\30p5ms';
name_pattern1 = "H2B_DMSO"; name_pattern2 = "_Cell01";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20231203_1016_epoch5_Ch1";
data_struct(iSamp).SampleName = "dense H2B";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20231203_1016_epoch5_Ch1_SR_pred_v2.csv";
FinalResults{iSamp}.locPerFrame = [];

iSamp = 3;
rootPath = 'F:\PROCESS-SPT\2024\20240419_U2OS_1x2x3x6xNLS_mPALM';
name_pattern1 = "1xNLS"; name_pattern2 = "_Cell01";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20231203_1016_epoch5_Ch1";
data_struct(iSamp).SampleName = "sparse 1xNLS";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20231203_1016_epoch5_Ch1_SR_pred_v2.csv";
FinalResults{iSamp}.locPerFrame = [];

iSamp = 4;
rootPath = 'F:\PROCESS-SPT\2024\20240303_U2OS_1xNLS_PA646_MPALM_dyeOnGlass';
name_pattern1 = "1xNLS"; name_pattern2 = "_Cell01";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20231203_1016_epoch5_Ch1";
data_struct(iSamp).SampleName = "dense 1xNLS";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20231203_1016_epoch5_Ch1_SR_pred_v2.csv";
FinalResults{iSamp}.locPerFrame = [];

%% Test if mEosEM works for mPALM: biased towards bound

clearvars data_struct FinalResults
% general parameters
FitParams.DoSingleFit = false;
FitParams.DoSinglePlot = false;
FitParams.DoMergeFit = false;
FitParams.DoMergePlot = true;

% load sample
iSamp = 1;
rootPath = 'E:\Submission\2023_MotionBlur\PA646MBX_H2B_1x6xNLS_beads\1xNLS\UNet_mask_MBX_20231203_1016_epoch5';
name_pattern = "30p5ms_20kframe";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Result_UNet_mask_MBX_20231203_1016_epoch5.mat";
data_struct(iSamp).SampleName = "1xNLS-HaloPA646";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20231203_1016_epoch5_SR_pred_v2.csv";%"SR_pred.csv";

% load sample
iSamp = 2;
rootPath = 'E:\PROCESS-SPT\2024\20240522_1xNLS_mEosEM_MPALMtest\';
name_pattern = "20240419_U2OS_1xNLS_mEosEM_40pert_30p5ms_10kframe_005_Cell01";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20231203_1016_epoch5_Ch1.mat";
data_struct(iSamp).SampleName = "1xNLS-mEosEM_条件5";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20231203_1016_epoch5_Ch1_SR_pred_v2.csv";%"SR_pred.csv";

% load sample
iSamp = 3;
rootPath = 'E:\PROCESS-SPT\2024\20240522_1xNLS_mEosEM_MPALMtest\';
name_pattern = "20240419_U2OS_1xNLS_mEosEM_40pert_30p5ms_10kframe_005_Cell01";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20240620_2035_epoch20_Ch1.mat";
data_struct(iSamp).SampleName = "1xNLS-mEosEM_new条件5";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20240620_2035_epoch20_Ch1_SR_pred_v3.csv";%"SR_pred.csv";
% 
% 
% % load sample
% iSamp = 2;
% rootPath = 'E:\PROCESS-SPT\2024\20240522_1xNLS_mEosEM_MPALMtest\';
% name_pattern = "20240419_U2OS_1xNLS_mEosEM_40pert_30p5ms_10kframe_001_Cell01";
% temp_file = dir(fullfile(rootPath));
% temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern);
% data_struct(iSamp).workspaces = temp_file(idx);
% data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20231203_1016_epoch5_Ch1.mat";
% data_struct(iSamp).SampleName = "1xNLS-mEosEM_条件1";
% data_struct(iSamp).SRFile = "UNet_mask_MBX_20231203_1016_epoch5_Ch1_SR_pred_v2.csv";%"SR_pred.csv";
% 
% % load sample
% iSamp = 3;
% rootPath = 'E:\PROCESS-SPT\2024\20240522_1xNLS_mEosEM_MPALMtest\';
% name_pattern = "20240419_U2OS_1xNLS_mEosEM_50pert_30p5ms_10kframe_002_Cell01";
% temp_file = dir(fullfile(rootPath));
% temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern);
% data_struct(iSamp).workspaces = temp_file(idx);
% data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20231203_1016_epoch5_Ch1.mat";
% data_struct(iSamp).SampleName = "1xNLS-mEosEM_条件2";
% data_struct(iSamp).SRFile = "UNet_mask_MBX_20231203_1016_epoch5_Ch1_SR_pred_v2.csv";%"SR_pred.csv";
% 
% % load sample
% iSamp = 4;
% rootPath = 'E:\PROCESS-SPT\2024\20240522_1xNLS_mEosEM_MPALMtest\';
% name_pattern = "20240419_U2OS_1xNLS_mEosEM_60pert_30p5ms_10kframe_003_Cell01";
% temp_file = dir(fullfile(rootPath));
% temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern);
% data_struct(iSamp).workspaces = temp_file(idx);
% data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20231203_1016_epoch5_Ch1.mat";
% data_struct(iSamp).SampleName = "1xNLS-mEosEM_条件3";
% data_struct(iSamp).SRFile = "UNet_mask_MBX_20231203_1016_epoch5_Ch1_SR_pred_v2.csv";%"SR_pred.csv";
% 
% % load sample
% iSamp = 5;
% rootPath = 'E:\PROCESS-SPT\2024\20240522_1xNLS_mEosEM_MPALMtest\';
% name_pattern = "20240419_U2OS_1xNLS_mEosEM_100pert_30p5ms_10kframe_004_Cell01";
% temp_file = dir(fullfile(rootPath));
% temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern);
% data_struct(iSamp).workspaces = temp_file(idx);
% data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20231203_1016_epoch5_Ch1.mat";
% data_struct(iSamp).SampleName = "1xNLS-mEosEM_条件4";
% data_struct(iSamp).SRFile = "UNet_mask_MBX_20231203_1016_epoch5_Ch1_SR_pred_v2.csv";%"SR_pred.csv";
% 
% % load sample
% iSamp = 6;
% rootPath = 'E:\PROCESS-SPT\2024\20240522_1xNLS_mEosEM_MPALMtest\';
% name_pattern = "20240419_U2OS_1xNLS_mEosEM_40pert_30p5ms_10kframe_005_Cell01";
% temp_file = dir(fullfile(rootPath));
% temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern);
% data_struct(iSamp).workspaces = temp_file(idx);
% data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20231203_1016_epoch5_Ch1.mat";
% data_struct(iSamp).SampleName = "1xNLS-mEosEM_条件5";
% data_struct(iSamp).SRFile = "UNet_mask_MBX_20231203_1016_epoch5_Ch1_SR_pred_v2.csv";%"SR_pred.csv";
% 
% % load sample
% iSamp = 7;
% rootPath = 'E:\PROCESS-SPT\2024\20240522_1xNLS_mEosEM_MPALMtest\';
% name_pattern = "20240419_U2OS_1xNLS_mEosEM_40pert_30p5ms_10kframe_006_Cell01";
% temp_file = dir(fullfile(rootPath));
% temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern);
% data_struct(iSamp).workspaces = temp_file(idx);
% data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20231203_1016_epoch5_Ch1.mat";
% data_struct(iSamp).SampleName = "1xNLS-mEosEM_条件6";
% data_struct(iSamp).SRFile = "UNet_mask_MBX_20231203_1016_epoch5_Ch1_SR_pred_v2.csv";%"SR_pred.csv";

%% 1x2x3x6xNLS and H2B
clearvars data_struct FinalResults
% general parameters
FitParams.DoSingleFit = true;
FitParams.DoSinglePlot = false;
FitParams.DoMergeFit = true;
FitParams.DoMergePlot = true;

% load sample
iSamp = 1;
rootPath = 'G:\PROCESS-SPT\2024\20240419_U2OS_1x2x3x6xNLS_mPALM';
name_pattern1 = "1xNLS"; name_pattern2 = "_Cell01";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20231203_1016_epoch5_Ch1";
data_struct(iSamp).SampleName = "1xNLS";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20231203_1016_epoch5_Ch1_SR_pred_v2.csv";

iSamp = 2;
rootPath = 'G:\PROCESS-SPT\2024\20240419_U2OS_1x2x3x6xNLS_mPALM';
name_pattern1 = "2xNLS"; name_pattern2 = "_Cell01";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20231203_1016_epoch5_Ch1";
data_struct(iSamp).SampleName = "2xNLS";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20231203_1016_epoch5_Ch1_SR_pred_v2.csv";

iSamp = 3;
rootPath = 'G:\PROCESS-SPT\2024\20240419_U2OS_1x2x3x6xNLS_mPALM';
name_pattern1 = "3xNLS"; name_pattern2 = "_Cell01";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20231203_1016_epoch5_Ch1";
data_struct(iSamp).SampleName = "3xNLS";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20231203_1016_epoch5_Ch1_SR_pred_v2.csv";

iSamp = 4;
rootPath = 'G:\PROCESS-SPT\2024\20240419_U2OS_1x2x3x6xNLS_mPALM';
name_pattern1 = "6xNLS"; name_pattern2 = "_Cell01";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20231203_1016_epoch5_Ch1";
data_struct(iSamp).SampleName = "6xNLS";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20231203_1016_epoch5_Ch1_SR_pred_v2.csv";

iSamp = 5;
rootPath = 'G:\PROCESS-SPT\2024\20240602_U2OS_FOXA2_MPALM';
name_pattern1 = "30p5ms_5kframe"; name_pattern2 = "_Cell01";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20231203_1016_epoch5_Ch1";
data_struct(iSamp).SampleName = "FOXA2";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20231203_1016_epoch5_Ch1_SR_pred_v3.csv";%"SR_pred.csv";

iSamp = 6;
rootPath = 'G:\Submission\2023_MotionBlur\PA646MBX_H2B_1x6xNLS_beads\H2B\UNet_mask_MBX_20231203_1016_epoch5';
name_pattern = "30p5ms_5kframe";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Result_UNet_mask_MBX_20231203_1016_epoch5.mat";
data_struct(iSamp).SampleName = "H2B";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20231203_1016_epoch5_SR_pred_v2.csv";

iSamp = 7;
rootPath = 'G:\Submission\2023_MotionBlur\PA646MBX_H2B_1x6xNLS_beads\beads';
name_pattern = "30p5ms_5kframe_01";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Result_UNet_mask_MBX_20230712_1000_epoch12.mat";
data_struct(iSamp).SampleName = "beads";
data_struct(iSamp).SRFile = "SRpred_singleBeads_MBX_20231220_110nmPix_rep2_epoch9.csv";

%% OPRM1
clearvars data_struct FinalResults
% general parameters
FitParams.DoSingleFit = false;
FitParams.DoSinglePlot = false;
FitParams.DoMergeFit = true;
FitParams.DoMergePlot = true;

% load sample
iSamp = 1;
rootPath = 'E:\PROCESS-SPT\2024\20240524_U2OS_OPRM1_MPALM\';
name_pattern1 = "DMSO"; name_pattern2 = "_Cell01";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20231203_1016_epoch5_Ch1";
data_struct(iSamp).SampleName = "OPRM1 DMSO";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20231203_1016_epoch5_Ch1_SR_pred_v2.csv";

iSamp = 2;
rootPath = 'E:\PROCESS-SPT\2024\20240524_U2OS_OPRM1_MPALM\';
name_pattern1 = "DAMGO10uM"; name_pattern2 = "_Cell01";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20231203_1016_epoch5_Ch1";
data_struct(iSamp).SampleName = "OPRM1 DAMGO10uM";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20231203_1016_epoch5_Ch1_SR_pred_v2.csv";

%% H2B FLA treatment
clearvars data_struct FinalResults
% general parameters
FitParams.DoSingleFit = true;
FitParams.DoSinglePlot = false;
FitParams.DoMergeFit = false;
FitParams.DoMergePlot = true;

% load sample
iSamp = 1;
rootPath = 'E:\PROCESS-SPT\2024\20240320_U2OS_EF1a_H2B-Halo_PA646_DRB_FLA\30p5ms\';
name_pattern1 = "DMSO"; name_pattern2 = "_Cell01";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20231203_1016_epoch5_Ch1";
data_struct(iSamp).SampleName = "H2B DMSO";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20231203_1016_epoch5_Ch1_SR_pred_v2.csv";

iSamp = 2;
rootPath = 'E:\PROCESS-SPT\2024\20240320_U2OS_EF1a_H2B-Halo_PA646_DRB_FLA\30p5ms\';
name_pattern1 = "10uMFLA"; name_pattern2 = "_Cell01";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20231203_1016_epoch5_Ch1";
data_struct(iSamp).SampleName = "H2B FLA10uM";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20231203_1016_epoch5_Ch1_SR_pred_v2.csv";

%% TDP43 Heatshock treatment
clearvars data_struct FinalResults
% general parameters
FitParams.DoSingleFit = true;
FitParams.DoSinglePlot = false;
FitParams.DoMergeFit = false;
FitParams.DoMergePlot = true;

% load sample
iSamp = 1;
rootPath = 'E:\PROCESS-SPT\2023\20230917_HeLa_TDP43-SNAP_TMRBlink_HeatShock_mPALM\';
name_pattern1 = "TDP43_100uMTMRB"; name_pattern2 = "_Cell01";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20231203_1016_epoch5_Ch1";
data_struct(iSamp).SampleName = "TDP43 Normal Temp";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20231203_1016_epoch5_Ch1_SR_pred_v2.csv";

iSamp = 2;
rootPath = 'E:\PROCESS-SPT\2023\20230917_HeLa_TDP43-SNAP_TMRBlink_HeatShock_mPALM\';
name_pattern1 = "TDP43AfterHS1hr_100uMTMRB"; name_pattern2 = "_Cell01";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20231203_1016_epoch5_Ch1";
data_struct(iSamp).SampleName = "TDP43 HS";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20231203_1016_epoch5_Ch1_SR_pred_v2.csv";

%% HSF1 Heatshock treatment
clearvars data_struct FinalResults
% general parameters
FitParams.DoSingleFit = true;
FitParams.DoSinglePlot = false;
FitParams.DoMergeFit = false;
FitParams.DoMergePlot = true;

% load sample
iSamp = 1;
rootPath = 'E:\PROCESS-SPT\WangBo\20231114_U2OS-HSF1-Halo-TDP43-SNAP_mPALM\';
name_pattern1 = "HSF1-Halo-TDP43-SNAP_642_25nMPA646_30p5ms_10kframe"; name_pattern2 = "_Cell01";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20231203_1016_epoch5_Ch1";
data_struct(iSamp).SampleName = "HSF1 Normal Temp";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20231203_1016_epoch5_Ch1_SR_pred_v2.csv";

iSamp = 2;
rootPath = 'E:\PROCESS-SPT\WangBo\20231114_U2OS-HSF1-Halo-TDP43-SNAP_mPALM\';
name_pattern1 = "HSF1-Halo-TDP43-SNAP_80minHS_642_25nMPA646_30p5ms_10kframe"; name_pattern2 = "_Cell01";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20231203_1016_epoch5_Ch1";
data_struct(iSamp).SampleName = "HSF1 HS";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20231203_1016_epoch5_Ch1_SR_pred_v2.csv";


%% 1xNLS/6xNLS/FOXA2/H2B/beads
clearvars data_struct FinalResults
% general parameters
FitParams.DoSingleFit = false;
FitParams.DoSinglePlot = false;
FitParams.DoMergeFit = false;
FitParams.DoMergePlot = true;
% hist_edges = linspace(0,4,100);

% load sample
iSamp = 1;
rootPath = 'E:\Submission\2023_MotionBlur\PA646MBX_H2B_1x6xNLS_beads\1xNLS\UNet_mask_MBX_20231203_1016_epoch5';
name_pattern = "30p5ms_20kframe";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Result_UNet_mask_MBX_20231203_1016_epoch5.mat";
data_struct(iSamp).SampleName = "1xNLS";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20231203_1016_epoch5_SR_pred_v2.csv";%"SR_pred.csv";

iSamp = 2;
rootPath = 'E:\Submission\2023_MotionBlur\PA646MBX_H2B_1x6xNLS_beads\6xNLS\UNet_mask_MBX_20231203_1016_epoch5';
name_pattern = "30p5ms_5kframe";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Result_UNet_mask_MBX_20231203_1016_epoch5.mat";
data_struct(iSamp).SampleName = "6xNLS";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20231203_1016_epoch5_SR_pred_v2.csv";%"SR_pred.csv";

iSamp = 3;
rootPath = 'E:\PROCESS-SPT\2024\20240602_U2OS_FOXA2_MPALM';
name_pattern1 = "30p5ms_5kframe"; name_pattern2 = "_Cell01";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20231203_1016_epoch5_Ch1";
data_struct(iSamp).SampleName = "FOXA2";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20231203_1016_epoch5_Ch1_SR_pred_v3.csv";%"SR_pred.csv";

iSamp = 4;
rootPath = 'E:\Submission\2023_MotionBlur\PA646MBX_H2B_1x6xNLS_beads\H2B\UNet_mask_MBX_20231203_1016_epoch5';
name_pattern = "30p5ms_5kframe";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Result_UNet_mask_MBX_20231203_1016_epoch5.mat";
data_struct(iSamp).SampleName = "H2B";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20231203_1016_epoch5_SR_pred_v2.csv";%"SR_pred.csv";

iSamp = 5;
rootPath = 'E:\Submission\2023_MotionBlur\PA646MBX_H2B_1x6xNLS_beads\beads';
name_pattern = "30p5ms_5kframe_01";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Result_UNet_mask_MBX_20230712_1000_epoch12.mat";
data_struct(iSamp).SampleName = "beads";
data_struct(iSamp).SRFile = "SRpred_singleBeads_MBX_20231220_110nmPix_rep2_epoch9.csv";

%% Compare UNet UNet_mask_MBX_20231203_1016_epoch5_Ch1 and UNet_mask_MBX_20240620_2035_epoch20_Ch1

clearvars data_struct FinalResults
% general parameters
FitParams.DoSingleFit = true;
FitParams.DoSinglePlot = false;
FitParams.DoMergeFit = true;
FitParams.DoMergePlot = true;

% load sample
iSamp = 1;
rootPath = 'E:\PROCESS-SPT\2024\20240419_U2OS_1x2x3x6xNLS_mPALM';
name_pattern1 = "1xNLS"; name_pattern2 = "_Cell01";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20231203_1016_epoch5_Ch1";
data_struct(iSamp).SampleName = "1xNLS_MBX_20231203_1016_epoch5_Ch1";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20231203_1016_epoch5_Ch1_SR_pred_v2.csv";

iSamp = 2;
rootPath = 'E:\PROCESS-SPT\2024\20240419_U2OS_1x2x3x6xNLS_mPALM';
name_pattern1 = "2xNLS"; name_pattern2 = "_Cell01";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20231203_1016_epoch5_Ch1";
data_struct(iSamp).SampleName = "2xNLS_MBX_20231203_1016_epoch5_Ch1";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20231203_1016_epoch5_Ch1_SR_pred_v2.csv";

iSamp = 3;
rootPath = 'E:\PROCESS-SPT\2024\20240419_U2OS_1x2x3x6xNLS_mPALM';
name_pattern1 = "3xNLS"; name_pattern2 = "_Cell01";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20231203_1016_epoch5_Ch1";
data_struct(iSamp).SampleName = "3xNLS_MBX_20231203_1016_epoch5_Ch1";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20231203_1016_epoch5_Ch1_SR_pred_v2.csv";

iSamp = 4;
rootPath = 'E:\PROCESS-SPT\2024\20240419_U2OS_1x2x3x6xNLS_mPALM';
name_pattern1 = "6xNLS"; name_pattern2 = "_Cell01";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20231203_1016_epoch5_Ch1";
data_struct(iSamp).SampleName = "6xNLS_MBX_20231203_1016_epoch5_Ch1";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20231203_1016_epoch5_Ch1_SR_pred_v2.csv";

iSamp = 5;
rootPath = 'E:\PROCESS-SPT\2024\20240419_U2OS_1x2x3x6xNLS_mPALM';
name_pattern1 = "1xNLS"; name_pattern2 = "_Cell01";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20240620_2035_epoch20_Ch1";
data_struct(iSamp).SampleName = "1xNLS_MBX_20240620_2035_epoch20_Ch1";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20240620_2035_epoch20_Ch1_SR_pred_v3.csv";

iSamp = 6;
rootPath = 'E:\PROCESS-SPT\2024\20240419_U2OS_1x2x3x6xNLS_mPALM';
name_pattern1 = "2xNLS"; name_pattern2 = "_Cell01";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20240620_2035_epoch20_Ch1";
data_struct(iSamp).SampleName = "2xNLS_MBX_20240620_2035_epoch20_Ch1";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20240620_2035_epoch20_Ch1_SR_pred_v3.csv";

iSamp = 7;
rootPath = 'E:\PROCESS-SPT\2024\20240419_U2OS_1x2x3x6xNLS_mPALM';
name_pattern1 = "3xNLS"; name_pattern2 = "_Cell01";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20240620_2035_epoch20_Ch1";
data_struct(iSamp).SampleName = "3xNLS_MBX_20240620_2035_epoch20_Ch1";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20240620_2035_epoch20_Ch1_SR_pred_v3.csv";

iSamp = 8;
rootPath = 'E:\PROCESS-SPT\2024\20240419_U2OS_1x2x3x6xNLS_mPALM';
name_pattern1 = "6xNLS"; name_pattern2 = "_Cell01";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20240620_2035_epoch20_Ch1";
data_struct(iSamp).SampleName = "6xNLS_MBX_20240620_2035_epoch20_Ch1";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20240620_2035_epoch20_Ch1_SR_pred_v3.csv";

%% Dual color RPB1 DMSO/THZ1

clearvars data_struct FinalResults
% general parameters
FitParams.DoSingleFit = true;
FitParams.DoSinglePlot = false;
FitParams.DoMergeFit = true;
FitParams.DoMergePlot = true;

% load sample
iSamp = 1;
rootPath = 'F:\PROCESS-SPT\2024\20240401_U2OS_RPB1Halo-H2BmEosEM_dualColor_MPALM';
name_pattern1 = "DMSO"; name_pattern2 = "_Cell01";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20240620_2035_epoch20_Ch2.mat";
data_struct(iSamp).SampleName = "RPB1Halo-H2BmEosEM DMSO";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20240620_2035_epoch20_Ch2_SR_pred_v3.csv";

iSamp = 2;
rootPath = 'F:\PROCESS-SPT\2024\20240401_U2OS_RPB1Halo-H2BmEosEM_dualColor_MPALM';
name_pattern1 = "1uMTHZ1"; name_pattern2 = "_Cell01";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20240620_2035_epoch20_Ch2.mat";
data_struct(iSamp).SampleName = "RPB1Halo-H2BmEosEM 1uM THZ1";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20240620_2035_epoch20_Ch2_SR_pred_v3.csv";

%% OPRM1

clearvars data_struct FinalResults
% general parameters
FitParams.DoSingleFit = false;
FitParams.DoSinglePlot = false;
FitParams.DoMergeFit = false;
FitParams.DoMergePlot = true;

% load sample
iSamp = 1;
rootPath = 'E:\PROCESS-SPT\2024\20240902_U2OS_TMsignal-OPRM1_PA646_MPALM';
name_pattern1 = "DMSO"; name_pattern2 = "_Cell01";name_pattern3 = "Clust01";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2) & contains(temp_Filenames,name_pattern3);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_NoLocFilter_UNet_mask_MBX_20240620_2035_epoch20_Ch1_Slice01.mat";
data_struct(iSamp).SampleName = "OPRM1-Halo DMSO 0902";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20240620_2035_epoch20_Ch1_SR_pred_v3.csv";

iSamp = 2;
rootPath = 'E:\PROCESS-SPT\2024\20240902_U2OS_TMsignal-OPRM1_PA646_MPALM';
name_pattern1 = "DAMGO1hr"; name_pattern2 = "_Cell01";name_pattern3 = "Clust01";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2) & contains(temp_Filenames,name_pattern3);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_NoLocFilter_UNet_mask_MBX_20240620_2035_epoch20_Ch1_Slice01.mat";
data_struct(iSamp).SampleName = "OPRM1-Halo 10uM DAMGO 1hr 0902";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20240620_2035_epoch20_Ch1_SR_pred_v3.csv";

iSamp = 3;
rootPath = 'E:\PROCESS-SPT\2024\20240902_U2OS_TMsignal-OPRM1_PA646_MPALM';
name_pattern1 = "NAX3hr"; name_pattern2 = "_Cell01";name_pattern3 = "Clust06";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2) & contains(temp_Filenames,name_pattern3);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_NoLocFilter_UNet_mask_MBX_20240620_2035_epoch20_Ch1_Slice01.mat";
data_struct(iSamp).SampleName = "OPRM1-Halo 10uM naxolone 3hr 0902";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20240620_2035_epoch20_Ch1_SR_pred_v3.csv";

iSamp = 4;
rootPath = 'E:\PROCESS-SPT\2024\20240904_U2OS_TMsignal-OPRM1_PA646_MPALM';
name_pattern1 = "DMSO"; name_pattern2 = "_Cell01";name_pattern3 = "Clust01";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2) & contains(temp_Filenames,name_pattern3);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20240620_2035_epoch20_Ch1_Slice01.mat";
data_struct(iSamp).SampleName = "OPRM1-Halo DMSO 0904";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20240620_2035_epoch20_Ch1_SR_pred_v3.csv";

iSamp = 5;
rootPath = 'E:\PROCESS-SPT\2024\20240904_U2OS_TMsignal-OPRM1_PA646_MPALM';
name_pattern1 = "1uMPZM21for1hr"; name_pattern2 = "_Cell01";name_pattern3 = "Clust02";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2) & contains(temp_Filenames,name_pattern3);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20240620_2035_epoch20_Ch1_Slice01.mat";
data_struct(iSamp).SampleName = "OPRM1-Halo 1uMPZM21 1hr 0904";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20240620_2035_epoch20_Ch1_SR_pred_v3.csv";

%% Mobility Spectrum using map

clearvars data_struct FinalResults
% general parameters
FitParams.DoSingleFit = false;
FitParams.DoSinglePlot = false;
FitParams.DoMergeFit = false;
FitParams.DoMergePlot = true;

% load sample
iSamp = 1;
rootPath = 'E:\PROCESS-SPT\2024\20240303_U2OS_1xNLS_PA646_MPALM_dyeOnGlass';
name_pattern1 = "1xNLS"; name_pattern2 = "_Cell01";name_pattern3 = "Clust02";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2) & contains(temp_Filenames,name_pattern3);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20240620_2035_epoch20_Ch1.mat";
data_struct(iSamp).SampleName = "1xNLS";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20240620_2035_epoch20_Ch1_SR_pred_v3.csv";

iSamp = 2;
rootPath = 'E:\PROCESS-SPT\2024\20240602_U2OS_FOXA2_MPALM';
name_pattern1 = "FOXA2"; name_pattern2 = "_Cell01";name_pattern3 = "Clust13";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2) & contains(temp_Filenames,name_pattern3);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20240620_2035_epoch20_Ch1.mat";
data_struct(iSamp).SampleName = "FOXA2";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20240620_2035_epoch20_Ch1_SR_pred_v3.csv";

iSamp = 3;
rootPath = 'E:\PROCESS-SPT\2023\20231013_U2OS_H2BOnly_HP1-H2B_mPALM\EF1a-H2B';
name_pattern1 = "H2B"; name_pattern2 = "_Cell01";name_pattern3 = "Clust01";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2) & contains(temp_Filenames,name_pattern3);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20240620_2035_epoch20_Ch1.mat";
data_struct(iSamp).SampleName = "H2B";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20240620_2035_epoch20_Ch1_SR_pred_v3.csv";

iSamp = 4;
rootPath = 'G:\PROCESS-SPT\20240712_U2OS_Paxillin_PA646_mPALM';
name_pattern1 = "Paxillin"; name_pattern2 = "_Cell01";name_pattern3 = "Clust01";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2) & contains(temp_Filenames,name_pattern3);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "  .mat";
data_struct(iSamp).SampleName = "paxillin";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20240620_epoch13_Ch1_SR_pred_v3.csv";

iSamp = 5;
rootPath = 'E:\PROCESS-SPT\2024\20240319_U2OS_NPM1_FBL_CD9';
name_pattern1 = "CD9"; name_pattern2 = "_Cell01";name_pattern3 = "Clust01";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2) & contains(temp_Filenames,name_pattern3);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20240620_2035_epoch20_Ch1.mat";
data_struct(iSamp).SampleName = "CD9";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20240620_2035_epoch20_Ch1_SR_pred_v3.csv";

%% Mobility Spectrum using map update using 1x2x3x6xNLS, FOXA2, H2B

clearvars data_struct FinalResults
% general parameters
FitParams.DoSingleFit = false;
FitParams.DoSinglePlot = false;
FitParams.DoMergeFit = false;
FitParams.DoMergePlot = true;

% load sample
iSamp = 1;
rootPath = 'F:\PROCESS-SPT\2024\20240303_U2OS_1xNLS_PA646_MPALM_dyeOnGlass';
name_pattern1 = "1xNLS"; name_pattern2 = "_Cell01";name_pattern3 = "Clust02";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2) & contains(temp_Filenames,name_pattern3);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20240620_2035_epoch20_Ch1.mat";
data_struct(iSamp).SampleName = "1xNLS";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20240620_2035_epoch20_Ch1_SR_pred_v3.csv";

iSamp = 2;
rootPath = 'F:\PROCESS-SPT\2024\20241015_U2OS_CMV_2x3x6xNLS-Halo_MPALM';
name_pattern1 = "2xNLS"; name_pattern2 = "_Cell01";name_pattern3 = "Clust03";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2) & contains(temp_Filenames,name_pattern3);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20240620_2035_epoch20_Ch1.mat";
data_struct(iSamp).SampleName = "2xNLS";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20240620_2035_epoch20_Ch1_SR_pred_v3.csv";

iSamp = 3;
rootPath = 'E:\PROCESS-SPT\2024\20241024_U2OS_CMV_2x3x6xNLS-Halo_MPALM';
name_pattern1 = "3xNLS"; name_pattern2 = "_Cell01";name_pattern3 = "Clust08";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2) & contains(temp_Filenames,name_pattern3);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20240620_2035_epoch20_Ch1.mat";
data_struct(iSamp).SampleName = "3xNLS";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20240620_2035_epoch20_Ch1_SR_pred_v3.csv";

iSamp = 4;
rootPath = 'E:\PROCESS-SPT\2024\20241024_U2OS_CMV_2x3x6xNLS-Halo_MPALM';
name_pattern1 = "6xNLS"; name_pattern2 = "_Cell01";name_pattern3 = "Clust13";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2) & contains(temp_Filenames,name_pattern3);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20240620_2035_epoch20_Ch1.mat";
data_struct(iSamp).SampleName = "6xNLS";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20240620_2035_epoch20_Ch1_SR_pred_v3.csv";

iSamp = 5;
rootPath = 'F:\PROCESS-SPT\2024\20240602_U2OS_FOXA2_MPALM';
name_pattern1 = "FOXA2"; name_pattern2 = "_Cell01";name_pattern3 = "Clust13";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2) & contains(temp_Filenames,name_pattern3);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20240620_2035_epoch20_Ch1.mat";
data_struct(iSamp).SampleName = "FOXA2";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20240620_2035_epoch20_Ch1_SR_pred_v3.csv";

iSamp = 6;
rootPath = 'F:\PROCESS-SPT\2023\20231013_U2OS_H2BOnly_HP1-H2B_mPALM\EF1a-H2B';
name_pattern1 = "H2B"; name_pattern2 = "_Cell01";name_pattern3 = "Clust02";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2) & contains(temp_Filenames,name_pattern3);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20240620_2035_epoch20_Ch1.mat";
data_struct(iSamp).SampleName = "H2B";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20240620_2035_epoch20_Ch1_SR_pred_v3.csv";

%% HXSiR OPRM1 Same cell DMSO->DAMGO->Naxolone

clearvars data_struct FinalResults
% general parameters
FitParams.DoSingleFit = false;
FitParams.DoSinglePlot = false;
FitParams.DoMergeFit = false;
FitParams.DoMergePlot = true;

% load sample
iSamp = 1;
rootPath = 'E:\PROCESS-SPT\2024\20240916_U2OS_OPRM1_HXSiR_testMPALM\';
name_pattern1 = "DMSO"; name_pattern2 = "60p5ms"; name_pattern3 = "10kframe"; name_pattern4 = "Cell01";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2) & contains(temp_Filenames,name_pattern3) & contains(temp_Filenames,name_pattern4);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20240620_2035_epoch20_Ch1.mat";
data_struct(iSamp).SampleName = "DMSO";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20240620_2035_epoch20_Ch1_SR_pred_v3.csv";

iSamp = 2;
rootPath = 'E:\PROCESS-SPT\2024\20240916_U2OS_OPRM1_HXSiR_testMPALM\';
name_pattern1 = "10uMDAMGO"; name_pattern2 = "60p5ms"; name_pattern3 = "10kframe";name_pattern4 = "Cell01";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2) & contains(temp_Filenames,name_pattern3)& contains(temp_Filenames,name_pattern4);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20240620_2035_epoch20_Ch1.mat";
data_struct(iSamp).SampleName = "10uMDAMGO";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20240620_2035_epoch20_Ch1_SR_pred_v3.csv";

iSamp = 3;
rootPath = 'E:\PROCESS-SPT\2024\20240916_U2OS_OPRM1_HXSiR_testMPALM\';
name_pattern1 = "10uMNXO"; name_pattern2 = "003"; name_pattern3 = "10kframe";name_pattern4 = "Cell01";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2) & contains(temp_Filenames,name_pattern3)& contains(temp_Filenames,name_pattern4);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20240620_2035_epoch20_Ch1.mat";
data_struct(iSamp).SampleName = "10uMNXO_10min";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20240620_2035_epoch20_Ch1_SR_pred_v3.csv";

iSamp = 4;
rootPath = 'E:\PROCESS-SPT\2024\20240916_U2OS_OPRM1_HXSiR_testMPALM\';
name_pattern1 = "10uMNXO"; name_pattern2 = "004"; name_pattern3 = "10kframe";name_pattern4 = "Cell01";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2) & contains(temp_Filenames,name_pattern3)& contains(temp_Filenames,name_pattern4);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20240620_2035_epoch20_Ch1.mat";
data_struct(iSamp).SampleName = "10uMNXO_20min";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20240620_2035_epoch20_Ch1_SR_pred_v3.csv";

iSamp = 5;
rootPath = 'E:\PROCESS-SPT\2024\20240916_U2OS_OPRM1_HXSiR_testMPALM\';
name_pattern1 = "10uMNXO"; name_pattern2 = "005"; name_pattern3 = "10kframe";name_pattern4 = "Cell01";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2) & contains(temp_Filenames,name_pattern3)& contains(temp_Filenames,name_pattern4);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20240620_2035_epoch20_Ch1.mat";
data_struct(iSamp).SampleName = "10uMNXO_30min";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20240620_2035_epoch20_Ch1_SR_pred_v3.csv";
%% Different exposure time

clearvars data_struct FinalResults
% general parameters
FitParams.DoSingleFit = false;
FitParams.DoSinglePlot = false;
FitParams.DoMergeFit = false;
FitParams.DoMergePlot = true;

% load sample
iSamp = 1;
rootPath = 'E:\PROCESS-SPT\2024\20240916_U2OS_OPRM1_HXSiR_testMPALM\';
name_pattern1 = "DMSO"; name_pattern2 = "30p5ms"; name_pattern3 = "Cell01";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2) & contains(temp_Filenames,name_pattern3);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20240620_2035_epoch20_Ch1.mat";
data_struct(iSamp).SampleName = "30p5ms";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20240620_2035_epoch20_Ch1_SR_pred_v3.csv";

% load sample
iSamp = 2;
rootPath = 'E:\PROCESS-SPT\2024\20240916_U2OS_OPRM1_HXSiR_testMPALM\';
name_pattern1 = "DMSO"; name_pattern2 = "40p5ms"; name_pattern3 = "Cell01";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2) & contains(temp_Filenames,name_pattern3);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20240620_2035_epoch20_Ch1.mat";
data_struct(iSamp).SampleName = "40p5ms";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20240620_2035_epoch20_Ch1_SR_pred_v3.csv";

% load sample
iSamp = 3;
rootPath = 'E:\PROCESS-SPT\2024\20240916_U2OS_OPRM1_HXSiR_testMPALM\';
name_pattern1 = "DMSO"; name_pattern2 = "50p5ms"; name_pattern3 = "Cell01";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2) & contains(temp_Filenames,name_pattern3);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20240620_2035_epoch20_Ch1.mat";
data_struct(iSamp).SampleName = "50p5ms";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20240620_2035_epoch20_Ch1_SR_pred_v3.csv";

% load sample
iSamp = 4;
rootPath = 'E:\PROCESS-SPT\2024\20240916_U2OS_OPRM1_HXSiR_testMPALM\';
name_pattern1 = "DMSO"; name_pattern2 = "60p5ms"; name_pattern3 = "Cell01"; name_pattern4 = "5kframe"; 
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2) & contains(temp_Filenames,name_pattern3) & contains(temp_Filenames,name_pattern4);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20240620_2035_epoch20_Ch1.mat";
data_struct(iSamp).SampleName = "60p5ms";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20240620_2035_epoch20_Ch1_SR_pred_v3.csv";

% load sample
iSamp = 5;
rootPath = 'E:\PROCESS-SPT\2024\20240916_U2OS_OPRM1_HXSiR_testMPALM\';
name_pattern1 = "DMSO"; name_pattern2 = "70p5ms"; name_pattern3 = "Cell01";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2) & contains(temp_Filenames,name_pattern3);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20240620_2035_epoch20_Ch1.mat";
data_struct(iSamp).SampleName = "70p5ms";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20240620_2035_epoch20_Ch1_SR_pred_v3.csv";

% load sample
iSamp = 6;
rootPath = 'E:\PROCESS-SPT\2024\20240916_U2OS_OPRM1_HXSiR_testMPALM\';
name_pattern1 = "DMSO"; name_pattern2 = "80p5ms"; name_pattern3 = "Cell01";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2) & contains(temp_Filenames,name_pattern3);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20240620_2035_epoch20_Ch1.mat";
data_struct(iSamp).SampleName = "80p5ms";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20240620_2035_epoch20_Ch1_SR_pred_v3.csv";

% load sample
iSamp = 7;
rootPath = 'E:\PROCESS-SPT\2024\20240916_U2OS_OPRM1_HXSiR_testMPALM\';
name_pattern1 = "DMSO"; name_pattern2 = "90p5ms"; name_pattern3 = "Cell01";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2) & contains(temp_Filenames,name_pattern3);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20240620_2035_epoch20_Ch1.mat";
data_struct(iSamp).SampleName = "90p5ms";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20240620_2035_epoch20_Ch1_SR_pred_v3.csv";

% load sample
iSamp = 8;
rootPath = 'E:\PROCESS-SPT\2024\20240916_U2OS_OPRM1_HXSiR_testMPALM\';
name_pattern1 = "DMSO"; name_pattern2 = "100p5ms"; name_pattern3 = "Cell01";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2) & contains(temp_Filenames,name_pattern3);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20240620_2035_epoch20_Ch1.mat";
data_struct(iSamp).SampleName = "100p5ms";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20240620_2035_epoch20_Ch1_SR_pred_v3.csv";

%% CJ9 BRD4 JQ1 treatment

clearvars data_struct FinalResults
% general parameters
FitParams.DoSingleFit = false;
FitParams.DoSinglePlot = false;
FitParams.DoMergeFit = false;
FitParams.DoMergePlot = true;

% load sample
iSamp = 1;
rootPath = 'E:\PROCESS-SPT\2024\20240925_CJ9_HaloBRD4_dyeTest_PA646_HMSiR_BD566b\';
name_pattern1 = "DMSO"; name_pattern2 = "Clust001"; name_pattern3 = "Cell01";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2) & contains(temp_Filenames,name_pattern3);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20240620_2035_epoch20_Ch1.mat";
data_struct(iSamp).SampleName = "DMSO";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20240620_2035_epoch20_Ch1_SR_pred_v3.csv";

% load sample
iSamp = 2;
rootPath = 'E:\PROCESS-SPT\2024\20240925_CJ9_HaloBRD4_dyeTest_PA646_HMSiR_BD566b\';
name_pattern1 = "1uMJQ1"; name_pattern2 = "Clust04"; name_pattern3 = "Cell01";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2) & contains(temp_Filenames,name_pattern3);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20240620_2035_epoch20_Ch1.mat";
data_struct(iSamp).SampleName = "1uMJQ1";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20240620_2035_epoch20_Ch1_SR_pred_v3.csv";

%% MOR many cell DMSO, DAMGO, Naxolone

clearvars data_struct FinalResults
% general parameters
FitParams.DoSingleFit = false;
FitParams.DoSinglePlot = false;
FitParams.DoMergeFit = false;
FitParams.DoMergePlot = true;

% load sample
iSamp = 1;
rootPath = 'E:\PROCESS-SPT\2024\20240902_U2OS_TMsignal-OPRM1_PA646_MPALM\';
name_pattern1 = "_DMSO_"; name_pattern2 = "Cell01";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20240620_2035_epoch20_Ch1.mat";
data_struct(iSamp).SampleName = "DMSO";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20240620_2035_epoch20_Ch1_SR_pred_v3.csv";

% load sample
iSamp = 2;
rootPath = 'E:\PROCESS-SPT\2024\20240902_U2OS_TMsignal-OPRM1_PA646_MPALM\';
name_pattern1 = "_DAMGO1uM_"; name_pattern2 = "Cell01";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2) & contains(temp_Filenames,name_pattern3);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20240620_2035_epoch20_Ch1.mat";
data_struct(iSamp).SampleName = "1uMJQ1";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20240620_2035_epoch20_Ch1_SR_pred_v3.csv";

%% 100ms data MOR single cell DMSO, DAMGO, Naxolone, PZM21 from Zuhui 

clearvars data_struct FinalResults
% general parameters
FitParams.DoSingleFit = false;
FitParams.DoSinglePlot = false;
FitParams.DoMergeFit = false;
FitParams.DoMergePlot = true;

% load sample
iSamp = 1;
rootPath = 'C:\Users\Zuhui\Downloads\20241006_U2OS_TMsignal-OPRM1_PA646_MPALM100\';
name_pattern1 = "DMSO"; name_pattern2 = "Cell01"; name_pattern3 = "Clust08";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2) & contains(temp_Filenames,name_pattern3);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20240620_2035_epoch20_Ch1.mat";
data_struct(iSamp).SampleName = "DMSO";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20240620_2035_epoch20_Ch1_SR_pred_v3.csv";

% load sample
iSamp = 2;
rootPath = 'C:\Users\Zuhui\Downloads\20241006_U2OS_TMsignal-OPRM1_PA646_MPALM100\';
name_pattern1 = "10uMDAMGO"; name_pattern2 = "Cell01"; name_pattern3 = "Clust04";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2)  & contains(temp_Filenames,name_pattern3);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20240620_2035_epoch20_Ch1.mat";
data_struct(iSamp).SampleName = "10uM DAMGO";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20240620_2035_epoch20_Ch1_SR_pred_v3.csv";

% load sample
iSamp = 3;
rootPath = 'C:\Users\Zuhui\Downloads\20241006_U2OS_TMsignal-OPRM1_PA646_MPALM100\';
name_pattern1 = "10uMNAX"; name_pattern2 = "Cell01"; name_pattern3 = "Clust07";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2) & contains(temp_Filenames,name_pattern3);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20240620_2035_epoch20_Ch1.mat";
data_struct(iSamp).SampleName = "10uM naxolone";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20240620_2035_epoch20_Ch1_SR_pred_v3.csv";

% load sample
iSamp = 4;
rootPath = 'C:\Users\Zuhui\Downloads\20241006_U2OS_TMsignal-OPRM1_PA646_MPALM100\';
name_pattern1 = "1uMPZM21"; name_pattern2 = "Cell01"; name_pattern3 = "Clust11";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2) & contains(temp_Filenames,name_pattern3);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20240620_2035_epoch20_Ch1.mat";
data_struct(iSamp).SampleName = "1uM PZM21";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20240620_2035_epoch20_Ch1_SR_pred_v3.csv";

%% 100ms data MOR many cell DMSO, DAMGO, Naxolone from Zuhui Wang

clearvars data_struct FinalResults
% general parameters
FitParams.DoSingleFit = false;
FitParams.DoSinglePlot = true;
FitParams.DoMergeFit = false;
FitParams.DoMergePlot = true;

% load sample
iSamp = 1;
rootPath = 'C:\Users\Zuhui\Downloads\20241006_U2OS_TMsignal-OPRM1_PA646_MPALM100\';
name_pattern1 = "DMSO"; name_pattern2 = "Cell01"; 
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2) ;
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20240620_2035_epoch20_Ch1.mat";
data_struct(iSamp).SampleName = "DMSO";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20240620_2035_epoch20_Ch1_SR_pred_v3.csv";

% load sample
iSamp = 2;
rootPath = 'C:\Users\Zuhui\Downloads\20241006_U2OS_TMsignal-OPRM1_PA646_MPALM100\';
name_pattern1 = "10uMDAMGO"; name_pattern2 = "Cell01"; 
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2)  ;
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20240620_2035_epoch20_Ch1.mat";
data_struct(iSamp).SampleName = "10uM DAMGO";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20240620_2035_epoch20_Ch1_SR_pred_v3.csv";

% load sample
iSamp = 3;
rootPath = 'C:\Users\Zuhui\Downloads\20241006_U2OS_TMsignal-OPRM1_PA646_MPALM100\';
name_pattern1 = "10uMNAX"; name_pattern2 = "Cell01"; 
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2) ;
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20240620_2035_epoch20_Ch1.mat";
data_struct(iSamp).SampleName = "10uM naxolone";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20240620_2035_epoch20_Ch1_SR_pred_v3.csv";

% load sample
iSamp = 4;
rootPath = 'C:\Users\Zuhui\Downloads\20241006_U2OS_TMsignal-OPRM1_PA646_MPALM100\';
name_pattern1 = "1uMPZM21"; name_pattern2 = "Cell01";
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2) ;
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20240620_2035_epoch20_Ch1.mat";
data_struct(iSamp).SampleName = "1uM PZM21";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20240620_2035_epoch20_Ch1_SR_pred_v3.csv";

%% 100ms data MOR many cell DMSO, DAMGO, Naxolone from LiuYiwen

clearvars data_struct FinalResults
% general parameters
FitParams.DoSingleFit = false;
FitParams.DoSinglePlot = true;
FitParams.DoMergeFit = false;
FitParams.DoMergePlot = true;

% load sample
iSamp = 1;
rootPath = 'G:\PROCESS-SPT\20241010_U2OS_OPRM1_PA646_MPALM\';
name_pattern1 = "DMSO"; name_pattern2 = "Cell01"; 
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20240620_epoch13_Ch1.mat";
data_struct(iSamp).SampleName = "DMSO";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20240620_epoch13_Ch1_SR_pred_v3.csv";

% load sample
iSamp = 2;
rootPath = 'G:\PROCESS-SPT\20241010_U2OS_OPRM1_PA646_MPALM\';
name_pattern1 = "DAMGO"; name_pattern2 = "Cell01"; 
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20240620_epoch13_Ch1.mat";
data_struct(iSamp).SampleName = "10uM DAMGO";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20240620_epoch13_Ch1_SR_pred_v3.csv";

% load sample
iSamp = 3;
rootPath = 'G:\PROCESS-SPT\20241010_U2OS_OPRM1_PA646_MPALM\';
name_pattern1 = "NAX"; name_pattern2 = "Cell01"; 
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20240620_epoch13_Ch1.mat";
data_struct(iSamp).SampleName = "10uM naxolone";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20240620_epoch13_Ch1_SR_pred_v3.csv";

% load sample
iSamp = 4;
rootPath = 'G:\PROCESS-SPT\20241010_U2OS_OPRM1_PA646_MPALM\';
name_pattern1 = "PZM21"; name_pattern2 = "Cell01"; 
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2);
data_struct(iSamp).workspaces = temp_file(idx);
data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20240620_epoch13_Ch1.mat";
data_struct(iSamp).SampleName = "1uM PZM21";
data_struct(iSamp).SRFile = "UNet_mask_MBX_20240620_epoch13_Ch1_SR_pred_v3.csv";


%% 100ms data MOR many cell DMSO, DAMGO, Naxolone from LiuYiwen and Zuhui merge

clearvars data_struct FinalResults
% general parameters
FitParams.DoSingleFit = true;
FitParams.DoSinglePlot = false;
FitParams.DoMergeFit = true;
FitParams.DoMergePlot = true;

% load sample
iSamp = 1;
rootPath = 'F:\PROCESS-SPT\2024\20241006-1010_U2OS_TMsignal-OPRM1_PA646_MPALM100\';
name_pattern1 = "DMSO"; name_pattern2 = "Cell01"; 
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2);
data_struct(iSamp).workspaces = temp_file(idx);
% data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20240620_epoch13_Ch1.mat";
data_struct(iSamp).SampleName = "DMSO";

% load sample
iSamp = 2;
rootPath = 'F:\PROCESS-SPT\2024\20241006-1010_U2OS_TMsignal-OPRM1_PA646_MPALM100\';
name_pattern1 = "DAMGO"; name_pattern2 = "Cell01"; 
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2);
data_struct(iSamp).workspaces = temp_file(idx);
% data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20240620_epoch13_Ch1.mat";
data_struct(iSamp).SampleName = "10uM DAMGO";

% load sample
iSamp = 3;
rootPath = 'F:\PROCESS-SPT\2024\20241006-1010_U2OS_TMsignal-OPRM1_PA646_MPALM100\';
name_pattern1 = "PZM21"; name_pattern2 = "Cell01"; 
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2);
data_struct(iSamp).workspaces = temp_file(idx);
% data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20240620_epoch13_Ch1.mat";
data_struct(iSamp).SampleName = "1uM PZM21";

% load sample
iSamp = 4;
rootPath = 'F:\PROCESS-SPT\2024\20241006-1010_U2OS_TMsignal-OPRM1_PA646_MPALM100\';
name_pattern1 = "NAX"; name_pattern2 = "Cell01"; 
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2);
data_struct(iSamp).workspaces = temp_file(idx);
% data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20240620_epoch13_Ch1.mat";
data_struct(iSamp).SampleName = "10uM naxolone";

%% H2B RDF and L(r)-r

clearvars data_struct FinalResults
% general parameters
FitParams.DoSingleFit = false;
FitParams.DoSinglePlot = false;
FitParams.DoMergeFit = false;
FitParams.DoMergePlot = false;

% load sample
iSamp = 1;
rootPath = 'G:\PROCESS-SPT\2023\20231013_U2OS_H2BOnly_HP1-H2B_mPALM\EF1a-H2B\';
name_pattern1 = "20kframe_02"; name_pattern2 = "Cell01"; 
temp_file = dir(fullfile(rootPath));
temp_Filenames = {temp_file.name}; idx = contains(temp_Filenames,name_pattern1) & contains(temp_Filenames,name_pattern2);
data_struct(iSamp).workspaces = temp_file(idx);
% data_struct(iSamp).UNetMat = "Blurdata_UNet_mask_MBX_20240620_epoch13_Ch1.mat";
data_struct(iSamp).SampleName = "H2B_Cell02";
