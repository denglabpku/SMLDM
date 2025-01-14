% IMPORTANT NOTE: raw_table.Xpos and raw_table.Ypos unit are nm below.

%% ================== 2023/08/17 newer version of mobility-PALM rendering without cluster analysis ===================== %%
close all;clc;clear;
% load the dataset
root_path = 'E:\PROCESS-SPT\2023\20231013_U2OS_H2BOnly_HP1-H2B_mPALM\EF1a-H2B\20231013_Clust02_EF1a-H2B_642_25uMPA646_30p5ms_20kframe_02_Cell01';
raw_table = readtable(fullfile(root_path,'UNet_mask_MBX_20231203_1016_epoch5_SR_pred_v2.csv'));
load(fullfile(root_path,'Blurdata_UNet_mask_MBX_20231203_1016_epoch5.mat'),'impars')

G = findgroups(raw_table.Frame);
locsPerFrame = splitapply(@numel, raw_table.PSF_idx,G);
fprintf("Average detections per frame %.1f \n",mean(locsPerFrame));
fprintf("Max detections per frame %.1f \n",max(locsPerFrame));
fprintf("Min detections per frame %.1f \n",min(locsPerFrame));

% MPALM_RENDERv3 render
RenderParam.camPixelSize_Ch1 = 110;
RenderParam.D_SRArea_const = 1162;
RenderParam.render_minlogD = -1;
RenderParam.render_maxlogD = 1;
RenderParam.dxy = 30; % original set 50
RenderParam.sigma_render = 30; % original set 50
RenderParam.DensityCutoff_densityPALM = [0 5000];
RenderParam.DensityCutoff_densityMPALM = [0 3000];
RenderParam.Xrange = [0 RenderParam.camPixelSize_Ch1*impars.ImWidth];
RenderParam.Yrange = [0 RenderParam.camPixelSize_Ch1*impars.ImHeight];

% load localization table
raw_table = raw_table(raw_table.TotalPhoton <= Inf,:); % filter out abnormal molecules
% raw_table = tab_psf_fitresult;
raw_table.Xpos = raw_table.Xpos * RenderParam.camPixelSize_Ch1;
raw_table.Ypos = raw_table.Ypos * RenderParam.camPixelSize_Ch1;

%% ========================= create a time series mobility-PALM movie ===================================== 

movie_savepath = 'C:\Users\Zuhui\Downloads\test';

% Create a VideoWriter object for a new video file. Use the 'Archival' profile to specify a Motion JPEG 2000 file with lossless compression.
v1 = VideoWriter(fullfile(movie_savepath,'densityPALM.avi'),'Motion JPEG AVI');
v1.Quality = 95;
v1.FrameRate = 5;

v2 = VideoWriter(fullfile(movie_savepath,'MPALM.avi'),'Motion JPEG AVI');
v2.Quality = 95;
v2.FrameRate = 5;

v3 = VideoWriter(fullfile(movie_savepath,'densityMPALM.avi'),'Motion JPEG AVI');
v3.Quality = 95;
v3.FrameRate = 5;

% Verify the type of video compression for the new file.
v1.VideoCompressionMethod

% Open the video file for writing. Then, write the image data in A to the file.
open(v1);open(v2);open(v3);
for iframe = 1:2000:20000 % 164*30.5/1000 every 5 sec
    trunc_raw_table = raw_table(raw_table.Frame >= iframe & raw_table.Frame < iframe+2000,:);
    [~,~,~,fig1,fig2,fig3] = MPALM_RENDERv3(trunc_raw_table,RenderParam);
    f_densityPALM = getframe(fig1);
    f_MPALM = getframe(fig2);
    f_densityMPALM = getframe(fig3);
    writeVideo(v1,f_densityPALM);
    writeVideo(v2,f_MPALM);
    writeVideo(v3,f_densityMPALM);
    close(fig1);
    close(fig2);
    close(fig3);
end

% Close the video file.
close(v1);
close(v2);
close(v3);




