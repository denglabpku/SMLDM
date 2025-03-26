% Script to simulate snapshot of molecules using simSPT trajectory and
% save noisy snapshot and binary mask molecule image pairs

%% DESCRIPTION
% 1) Simulate snapshot of molecules using simSPT trajectory (simulation
% time gap 1ms). One snapshot is generated from one trajectory with
% desired track length from simSPT. Then, save the snapshot and sharp
% (immobile) molecule image pair into one image and can be used as the
% ground truth datasets for training networks. The pixelized image pair
% have the exactly 0.11um/pix size.
% 2) A Gaussian noise and Possion modeled shot-noise is added into simulated snapshot.

% Obtain your microscopy PSF's standard deviation from extract_psf.ipynb.
% Dependency: bfmatlab >6.6.0; msdanalyzer from matlab add-on market.

% Reference:
% simSPT: "https://gitlab.com/tjiandarzacq-lab/simSPT" from Hansen, Anders S., Iryna Pustova, Claudia Cattoglio, Robert Tjian, and Xavier Darzacq. “CTCF and Cohesin Regulate Chromatin Loop Stability with Distinct Dynamics.” ELife 6 (May 3, 2017): e25776. https://doi.org/10.7554/eLife.25776. 
% msdanalyzer: Jean-Yves Tinevez (2025). Mean square displacement analysis of particles trajectories (https://github.com/tinevez/msdanalyzer), GitHub. Retrieved January 11, 2025.

%% Zuhui Wang
%% 2025/01/10
%% INTENTIONALLY KEEP BLANK
%% INTENTIONALLY KEEP BLANK
%% INTENTIONALLY KEEP BLANK

close all; clear;clc;
addpath(genpath('./SMLMM/step0_simSnapshot'),'-end')
addpath(genpath('./MATLAB/msdanalyzer'),'-end')

% change your desired imaging parameters
target_pixelSize = 0.11; % desired simulated camera pixel size 0.11 or 0.16, um
background_mean = [1000 3000]; %1152; %500; %1761;% 1152;
background_std = [400 600];%283; %40; %409;% 283;
good_snr_range = [19 35]; 
save_MSD = false; % save MSD of each trajectories
PhotonConversionRate = 0.9398923112; % (0.5e/ADU, EMGain300, QE at 660nm:0.9398923112) (0.5e/ADU, EMGain300, QE at 582nm:0.9696909185)

% import fitting results from averaged beads
% Please follow the extract_psf.ipynb to obtain your PSF sigma value
beads_PixelSize=0.11; % um per pixel, pixel size of beads image
mean_r0 = 1.4141*beads_PixelSize; % replace 1.4141 with your microscope's sigma value


%% Simulate snapshot from simSPT
% import simSPT with time delay 1ms
% simSPT code: ./simSPT -D1=0.03 -D2=2.0 -p1=0.5 -p2=0.5 -sigma=0.031 -dt=0.001 -n_traj=100000 -file=../simPSF/20230309_D2_p5_dt1ms.csv -seed=0
simSPT_path = '/path/to/simSPT_generated_tracks';
simSPT_file = [...
    % >>>>>>>> No localization error >>>>>>>
    "20231202_D0p01_pure_10kTraj_dt1ms.csv",...
    "20231202_D0p1_pure_10kTraj_dt1ms.csv",...
    "20231202_D01_pure_10kTraj_dt1ms.csv",...
    "20231202_D02_pure_10kTraj_dt1ms.csv",...
    "20231202_D03_pure_10kTraj_dt1ms.csv",...
    "20231202_D04_pure_10kTraj_dt1ms.csv",...
    "20231202_D05_pure_10kTraj_dt1ms.csv",...
    "20231202_D06_pure_30kTraj_dt1ms.csv",...
    "20231202_D07_pure_30kTraj_dt1ms.csv",...
    "20231202_D08_pure_30kTraj_dt1ms.csv",...
    "20231202_D09_pure_30kTraj_dt1ms.csv",...
    "20231202_D10_pure_30kTraj_dt1ms.csv",...
    % "20231202_logDm3p0_pure_10kTraj_dt1ms.csv",...
    % "20231202_logDm2p8_pure_10kTraj_dt1ms.csv",...
    % "20231202_logDm2p6_pure_10kTraj_dt1ms.csv",...
    % "20231202_logDm2p4_pure_10kTraj_dt1ms.csv",...
    % "20231202_logDm2p2_pure_10kTraj_dt1ms.csv",...
    % "20231202_logDm2p0_pure_10kTraj_dt1ms.csv",...
    % "20231202_logDm1p8_pure_10kTraj_dt1ms.csv",...
    % "20231202_logDm1p6_pure_10kTraj_dt1ms.csv",...
    % "20231202_logDm1p4_pure_10kTraj_dt1ms.csv",...
    % "20231202_logDm1p2_pure_10kTraj_dt1ms.csv",...
    % "20231202_logDm1p0_pure_10kTraj_dt1ms.csv",...
    % "20231202_logDm0p8_pure_10kTraj_dt1ms.csv",...
    % "20231202_logDm0p6_pure_10kTraj_dt1ms.csv",...
    % "20231202_logDm0p4_pure_10kTraj_dt1ms.csv",...
    % "20231202_logDm0p2_pure_10kTraj_dt1ms.csv",...
    % "20231202_logD0p0_pure_10kTraj_dt1ms.csv",...
    % "20231202_logD0p2_pure_10kTraj_dt1ms.csv",...
    % "20231202_logD0p4_pure_10kTraj_dt1ms.csv",...
    % "20231202_logD0p6_pure_10kTraj_dt1ms.csv",...
    % "20231202_logD0p8_pure_30kTraj_dt1ms.csv",...
    % "20231202_logD0p9_pure_30kTraj_dt1ms.csv",...
    % "20231202_logD1p0_pure_30kTraj_dt1ms.csv",...
    % "20231202_logD1p1_pure_30kTraj_dt1ms.csv",...
    % "20231202_logD1p2_pure_30kTraj_dt1ms.csv",...
    % "20231202_logD1p3_pure_50kTraj_dt1ms.csv",...
    % "20231202_logD1p4_pure_50kTraj_dt1ms.csv",...
    % "20231202_logD1p5_pure_50kTraj_dt1ms.csv",...
    % <<<<<<<<< No localization error <<<<<<<<<<
    ];

simSPT_D = [...
    "p01", "p1","01","02","03","04","05","06","07","08","09", "10",...
    % "logDm3p0",...
    % "logDm2p8",...
    % "logDm2p6",...
    % "logDm2p4",...
    % "logDm2p2",...
    % "logDm2p0",...
    % "logDm1p8",...
    % "logDm1p6",...
    % "logDm1p4",...
    % "logDm1p2",...
    % "logDm1p0",...
    % "logDm0p8",...
    % "logDm0p6",...
    % "logDm0p4",...
    % "logDm0p2",...
    % "logD0p0",...
    % "logD0p2",...
    % "logD0p4",...
    % "logD0p6",...
    % "logD0p8",...
    % "logD0p9",...
    % "logD1p0",...
    % "logD1p1",...
    % "logD1p2",...
    % "logD1p3",...
    % "logD1p4",...
    % "logD1p5",...
    ];

% DiffCoeff = [0.00100    0.00158	0.00251	0.00398	0.00631	0.01000	0.01585	0.02512	0.03981	0.06310	0.10000	0.15849	0.25119	0.39811	0.63096	1.00000	1.58489	2.51189	3.98107	6.30957	7.94328	10.00000	12.58925	15.84893	19.95262	25.11886	31.62278];
DiffCoeff = [0.01 0.1 1:10];

% Starting index of simSPT trajectories
StartIdx = 1;
% maximal number of generated image pair at given D and exposure time
max_imgpair = 110; 
% desired exposure time to generate PSF
target_exposure = 30; %[1:10 20:10:90]; unit ms

% parent directory to save image pairs
img_pair_saveDir = '/path/to/save/simulated_dataset_forUNet';
mkdir(fullfile(img_pair_saveDir));
mkdir(fullfile(img_pair_saveDir,'imgs'));
mkdir(fullfile(img_pair_saveDir,'masks'));
% mkdir(fullfile(img_pair_saveDir,'tracks'));
% mkdir(fullfile(img_pair_saveDir,'trackHeatmap'));
% mkdir(fullfile(img_pair_saveDir,'locHeatmap'));

% save MSD of each images
if save_MSD
    MSD_savePath = img_pair_saveDir; % path to save MSD results
    MSD_T = table('Size',[0,9],'VariableTypes',{'cell','double','double','double','double','double','double','double','double'},...
        'VariableNames',{'ImageName','DiffCoeff','ExposureTime','SNR','background','noise','new_track_id','MSD','rep_id'});
    writetable(MSD_T,fullfile(MSD_savePath,'sim_summary.csv'),'WriteRowNames',true);
end

% initiate the log file to save image pair info
fileID = fopen(fullfile(img_pair_saveDir,'Summary.txt'),'w');
fprintf(fileID,'%s\t%g\n','target_pixelSize',target_pixelSize);
fprintf(fileID,'%s\t%d\t%d\n','background_mean',background_mean(1),background_mean(2));
fprintf(fileID,'%s\t%d\t%d\n','background_std',background_std(1),background_std(2));
fprintf(fileID,'%s\t%d\t%d\n','good_snr_range',good_snr_range(1),good_snr_range(2));
fprintf(fileID,'%s\t%s\t%s\n','simSPT_D','exposureTime','imagePair#');

%% ================ Start snapshot simulation ================ %%

for iSPT = 1:length(simSPT_file)
    fprintf('Processing simSPT_file: %s\n',simSPT_file{iSPT})
    t_raw = readtable(fullfile(simSPT_path,simSPT_file{iSPT}));
    
    % convert xy coord into pixel
    t_raw.xPixel = t_raw.x/target_pixelSize; % 0.11 um/pixel
    t_raw.yPixel = t_raw.y/target_pixelSize;
    
    % check track length
    [track_length, track_id, ~] = groupcounts(t_raw.trajectory);
        
    for idT = 1:length(target_exposure)
        fprintf('Processing target exposure time: %d ms ,simSPT file %d\n',target_exposure(idT),iSPT);
        
        % keep desired track length to have exactly target_exposure length. If track is long, truncate it using sliding window
        windowSize = 5;
        long_t = table();
        stopAllLoops = false; % Initialize flag variable
        total_subtrack_idx = 1;
        for i = 1:length(track_id)
            trajIndex = track_id(i);
            trajData = t_raw(t_raw.trajectory == trajIndex, :);
            dataHeight = height(trajData);

            if dataHeight >= target_exposure(idT)
                for j = 1:windowSize:dataHeight
                    if j+target_exposure(idT)-1 < dataHeight % not exceed the data height
                        subTrajData = trajData(j:j+target_exposure(idT)-1, :);
                        subTrajData.trajecotry_newIdx = ones(target_exposure(idT),1)*total_subtrack_idx;
                        long_t = vertcat(long_t,subTrajData);
                        total_subtrack_idx = total_subtrack_idx +1;
                        if total_subtrack_idx > max_imgpair
                            stopAllLoops = true; % Set flag variable to true
                            break; % Break out of innermost loop only
                        end
                    end
                end
                if stopAllLoops % Check flag variable
                    break; % If flag variable is true, break out of all nested loops
                end
            end
        end

        [long_track_length, long_track_id, ~] = groupcounts(long_t.trajecotry_newIdx);
        
        %% merge localization to generate psf in one pseudo-frame
        fprintf('Merging localization to generate psf ... \n');
        
        now_track_cell = {};
        % progressbarText(0);
        for iter = 1:length(long_track_id)
            % progressbarText(iter/length(long_track_id));
            now_track_cell{iter} = long_t(long_t.trajecotry_newIdx == long_track_id(iter),:);
        end

        if length(now_track_cell) < max_imgpair
            fprintf('Current trajectory# in this file is lower than max_imgpair!');
        end
        
        % initiate variable
        pix_filename_good = {};
        pix_filename_bad = {};

        good_noisy_pix_merge_z_snr = nan(length(now_track_cell),1);
        good_noisy_pix_merge_z_bkmean = nan(length(now_track_cell),1);
        good_noisy_pix_merge_z_bksigma = nan(length(now_track_cell),1);
        % bad_noisy_pix_merge_z_snr = nan(length(now_track_cell),1);
        % bad_noisy_pix_merge_z_bkmean = nan(length(now_track_cell),1);
        % bad_noisy_pix_merge_z_bksigma = nan(length(now_track_cell),1);
        % sum_jump_square = nan(length(now_track_cell),1);

        % parfor_progress(length(now_track_cell));
        progressbarText(0);
        for iter = 1:100
            progressbarText(iter/length(now_track_cell));
            now_track_cell{1} = now_track_cell{iter};
            now_track_cell{2} = now_track_cell{iter+1};
            now_track_cell{3} = now_track_cell{iter+2};
            now_track_cell{4} = now_track_cell{iter+3};

            rescale_pix_merge_z = zeros(32,32,4);
            mask_pix_merge_z = zeros(32,32,4);

            for itrack = 1:4
                now_track = now_track_cell{itrack};
                % simulate a snapshot by merge individual fitted psf from the same
                % trajectory
    
                % Calculate the range of x and y to span 3.52um (3.52/0.11=32
                % pixels) so that later can bin into 32x32 pixel with pixel
                % size 0.11um/px
    
                % Calculate the range of x and y to span 5.12um (5.12/0.16=32
                % pixels) so that later can bin into 32x32 pixel with pixel
                % size 0.16um/px
    
                range = 32*target_pixelSize/2;
                interval = 0.01; % spacing between x and y values
    
                if itrack == 1
                    x_min = mean(now_track.x) - 0.5*range;
                    x_max = mean(now_track.x) + 1.5*range;
                    y_min = mean(now_track.y) - 0.5*range;
                    y_max = mean(now_track.y) + 1.5*range;
                elseif itrack == 2
                    x_min = mean(now_track.x) - 1.5*range;
                    x_max = mean(now_track.x) + 0.5*range;
                    y_min = mean(now_track.y) - 0.5*range;
                    y_max = mean(now_track.y) + 1.5*range;
                elseif itrack == 3
                    x_min = mean(now_track.x) - 0.5*range;
                    x_max = mean(now_track.x) + 1.5*range;
                    y_min = mean(now_track.y) - 1.5*range;
                    y_max = mean(now_track.y) + 0.5*range;
                elseif itrack == 4
                    x_min = mean(now_track.x) - 1.5*range;
                    x_max = mean(now_track.x) + 0.5*range;
                    y_min = mean(now_track.y) - 1.5*range;
                    y_max = mean(now_track.y) + 0.5*range;
                end
                            
                % Generate x and y vectors that span 2xrange with 0.01 um gap size
                x = x_min:interval:x_max;
                y = y_min:interval:y_max;
    
                [Xgrid,Ygrid] = meshgrid(x,y);
            
                f = @(I,r0,m,x0,y0,x,y) I * exp(-(1/(2*r0^2))*((y - y0).^2 + (x - x0).^2)) + m;
            
                merge_z = zeros(length(Ygrid),length(Xgrid));
                % tempSharp_z = zeros(length(Ygrid),length(Xgrid),height(now_track)); % store instant frame intensity
    
                for i = 1:height(now_track)
                    
                    % noise_z = zeros(length(Ygrid),length(Xgrid)); % no noise consider the ideal situation
                    % rand_I = mean_I;
                    rand_I = 1;
                    % rand_I = normrnd(mean_I/10,sigma_I/100);
                    % rand_r0 = normrnd(mean_r0,sigma_r0);
                    rand_r0 = mean_r0;
                    rand_m = 0; % no background consider the ideal situation
            
                    % while (rand_I < 0 || rand_r0  < 0 || rand_m < 0 )
                    %     rand_I = normrnd(mean_I/10,sigma_I/100);
                    %     % rand_r0 = normrnd(mean_r0,sigma_r0);
                    %     rand_r0 = mean_r0;
                    %     rand_m = 0; 
                    % end
            
                    z = feval(f,rand_I,rand_r0,rand_m,now_track.x(i),now_track.y(i),Xgrid,Ygrid);
                    
                    merge_z = merge_z + z;
            
                end
    
                merge_z = merge_z/sum(merge_z(:));  % Normalize the matrix by dividing by the total volume
                
                % Grid the image to match the real camera pixel
                % Compute 32x32 bins with 0.11um bin width, i.e., pixel size is 0.11um/pix
                bin_width = target_pixelSize;%0.16;%0.11;
                num_bins = round(range * 2 / bin_width);
                Xedges = linspace(x_min, x_max, num_bins + 1);
                Yedges = linspace(y_min, y_max, num_bins + 1);
                
                % Bin Xgrid and Ygrid
                [N, pix_Xedge, pix_Yedge, binX, binY] = histcounts2(Xgrid, Ygrid, Xedges, Yedges);
                
                pix_merge_z = zeros(32,32);
                % Add the elements of C with the indices specified by A and B
                
                for i = 1:size(binX,1)
                    for j = 1:size(binX,2) 
                        pix_merge_z(binY(i,j),binX(i,j)) = merge_z(i,j) + pix_merge_z(binY(i,j),binX(i,j));
                    end
                end
    
                pix_merge_z = pix_merge_z/sum(pix_merge_z(:));  % Normalize the matrix by dividing by the total volume
    
                %========= save pixelized track image ======== %
                trackImg_pixelSize = 0.011; % unit um; pixel size of finer trajectory image
                trackImg_range = 320*trackImg_pixelSize/2;
                x_min = mean(now_track.x) - trackImg_range;
                x_max = mean(now_track.x) + trackImg_range;
                y_min = mean(now_track.y) - trackImg_range;
                y_max = mean(now_track.y) + trackImg_range;
                
                pix_track_x = round((now_track.x - x_min) / trackImg_pixelSize + 0.5);
                pix_track_y = round((now_track.y - y_min) / trackImg_pixelSize + 0.5);
    
                if min(pix_track_x) < 0 || min(pix_track_y) < 0 || max(pix_track_x) > 320 || max(pix_track_y) > 320 
    
                    continue % Current track cannot fit 320x320 image
    
                else
                                
                    %=========    create ground truth of motion blur/mask with good SNR ===============%
                    % % Create binary mask of motion blur, maxentropie only takes uint8 image
                    % Sort the normalized volume values in descending order
                    sorted_values = sort(pix_merge_z(:), 'descend');
    
                    % Find the index where the cumulative sum exceeds 95% of the total volume
                    cumulative_sum = cumsum(sorted_values);
                    threshold_index = find(cumulative_sum >= 0.95, 1);
    
                    % Create a binary mask using the threshold index
                    temp_mask_pix_merge_z = pix_merge_z >= sorted_values(threshold_index);
                    mask_pix_merge_z(:,:,itrack) = uint8(temp_mask_pix_merge_z); % convert from double to uint8
    
                    rescale_pix_merge_z(:,:,itrack) = rescale(pix_merge_z);
                end
            end
            rescale_pix_merge_z = sum(rescale_pix_merge_z,3);
            mask_pix_merge_z = sum(mask_pix_merge_z,3);
            mask_pix_merge_z = uint8(mask_pix_merge_z ~= 0);

            psf_num = bwlabel(mask_pix_merge_z,4);
            if max(psf_num,[],'all') < 4
                continue;
                % mask_pix_merge_z = uint8(zeros(size(mask_pix_merge_z)));
            end
                
            repIter = 1;
            while repIter <= 10
                % Adding noise to the pix_img_pair to mimic the real data
                % Create a 32x32 grayscale image of zeros
                % Set the range of the uniform distribution of SNR
                % SNR_range = [good_snr_range(1):good_snr_range(2)]; %[17 30]; 

                % Generate a single random number from the uniform distribution
                SNR = good_snr_range(1) + (good_snr_range(2)-good_snr_range(1))*rand(1,1);                
                bk_m = background_mean(1) + (background_mean(2)-background_mean(1))*rand(1,1);% background_mean;%1761;% 1152;
                bk_sig = background_std(1) + (background_std(2)-background_std(1))*rand(1,1);%background_std;%409;% 283;               

                % % % Noise model : Gaussian white noise
                % noisy_pix_merge_z = uint16(rescale_pix_merge_z*10^((SNR+20*log10(bk_sig))/20) + bk_m + bk_sig*randn(32,32));

                % Noise model : Poisson shot-noise + Gaussian white noise
                noisy_pix_merge_z = uint16(poissrnd(rescale_pix_merge_z*10^((SNR+20*log10(bk_sig))/20)) + bk_m + bk_sig*randn(32,32));

                % noise_img = uint16(bk_m + bk_sig*randn(32,32));
                good_noisy_pix_merge_z_snr(iter) = 10*log10((double(max(noisy_pix_merge_z,[],'all'))-double(median(noisy_pix_merge_z,'all')))^2/bk_sig^2);
                good_noisy_pix_merge_z_bkmean(iter) = bk_m;
                good_noisy_pix_merge_z_bksigma(iter) = bk_sig;

                % save the image mask and noisy snapshot
                imwrite(noisy_pix_merge_z,fullfile(img_pair_saveDir,'imgs',sprintf('0620_%s_TrackID%03d_D%s_SNR%d_dT%02dms_good_rep%03d.tif','Pixelized',now_track.trajecotry_newIdx(1),simSPT_D(iSPT),round(SNR),target_exposure(idT),repIter)));
                imwrite(mask_pix_merge_z,fullfile(img_pair_saveDir,'masks',sprintf('0620_%s_TrackID%03d_D%s_SNR%d_dT%02dms_good_rep%03d_mask.gif','Pixelized',now_track.trajecotry_newIdx(1),simSPT_D(iSPT),round(SNR),target_exposure(idT),repIter)));
                
                % % File Name that record SNR in table later
                % pix_filename_good{iter,1} = sprintf('%s_TrackID%012d_D%s_SNR%d_dT%02dms_good_rep%03d','Pixelized',now_track.trajecotry_newIdx(1),simSPT_D(iSPT),round(SNR),target_exposure(idT),repIter);
                %=========    create ground truth of motion blur/mask with good SNR ===============%
            
                % %========= save pixelized track image ======== %   
                % pix_track_img = zeros(320, 320);

                % % Draw lines connecting the pixels
                % for i = 2:length(pix_track_x)
                %     [x_line, y_line] = bresenham(pix_track_x(i-1), pix_track_y(i-1), pix_track_x(i), pix_track_y(i));
                %     ind = sub2ind([320, 320], y_line, x_line);
                %     pix_track_img(ind) = 1; % Set locations of the line to 1 in the pixelized image
                % end
                
                % % % (Debug only) Display the pixelized image
                % % figure;
                % % imagesc(pix_track_img);
                % % axis image ij;
                % % colormap(gray);
                % % imwrite(pix_track_img,fullfile(img_pair_saveDir,'tracks',sprintf('%s_TrackID%012d_D%s_SNR%d_dT%02dms_good_track.tif','Pixelized',now_track.trajecotry_newIdx(1),simSPT_D(iSPT),round(SNR),target_exposure(idT))));

                % % gaussian kernel standard deviation [pixels]
                % gaussian_sigma = 1;

                % % heatmap psf
                % psfHeatmap = fspecial('gauss',[7 7],gaussian_sigma);

                % % get the labels per frame in spikes and heatmaps
                % HeatmapImage = conv2(pix_track_img,psfHeatmap,'same');
                % imwrite(HeatmapImage,fullfile(img_pair_saveDir,'trackHeatmap',sprintf('%s_TrackID%012d_D%s_SNR%d_dT%02dms_good_trackHeatmap.tif','Pixelized',now_track.trajecotry_newIdx(1),simSPT_D(iSPT),round(SNR),target_exposure(idT))));
                
                % % % render localization heatmap
                % % pix_loc_img = zeros(320,320);
                % % pix_idx = sub2ind([320,320],pix_track_y,pix_track_x);
                % % pix_loc_img(pix_idx) = 1;
                % % render_pix_loc_img = conv2(pix_loc_img,psfHeatmap,'same');
                % % imwrite(render_pix_loc_img,fullfile(img_pair_saveDir,'locHeatmap',sprintf('%s_TrackID%012d_D%s_SNR%d_dT%02dms_good_locHeatmap.tif','Pixelized',now_track.trajecotry_newIdx(1),simSPT_D(iSPT),round(SNR),target_exposure(idT))));

                % %========= save pixelized track image ======== %
                
                %% calculate MSD of simSPT trajectories
                if save_MSD
                    % Initialize the msdanalyzer object
                    temp_simSPT_track_msd = msdanalyzer(2,'um','seconds');                      
                    temp_simSPT_track_msd = temp_simSPT_track_msd.addAll({[now_track.t,now_track.x,now_track.y]});
                    
                    % calculate the MSD of all imported tracks
                    temp_simSPT_track_msd = temp_simSPT_track_msd.computeMSD([]);
                    
                    % Append mean MSD into the table mean_T
                    ImageName = {sprintf('%s_TrackID%012d_D%s_SNR%d_dT%02dms_good_rep%03d.tif','Pixelized',now_track.trajecotry_newIdx(1),simSPT_D(iSPT),round(SNR),target_exposure(idT),repIter)};
                    new_line = table(ImageName,DiffCoeff(iSPT),target_exposure(idT),...
                        good_noisy_pix_merge_z_snr(iter), good_noisy_pix_merge_z_bkmean(iter), good_noisy_pix_merge_z_bksigma(iter),...
                        now_track.trajecotry_newIdx(1),temp_simSPT_track_msd.msd{1}(2,2),repIter); %MSD at 1*dT
                    writetable(new_line,fullfile(MSD_savePath,'sim_summary.csv'),'WriteMode','Append',...
                        'WriteVariableNames',false,'WriteRowNames',true); 
                end
                
                repIter = repIter+1;
            end
        end
        
        
        % save log 
        fprintf(fileID,'%3s\t%3d\t%12d\n',simSPT_D(iSPT),target_exposure(idT),length(now_track_cell));            
    end   
end

fclose(fileID);

% copy this script to the data folder as a log
sourceFile = mfilename('fullpath'); % Define the source and target file paths
copyfile([sourceFile '.m'], img_pair_saveDir);

disp("Simulated ground truth export finished!");
