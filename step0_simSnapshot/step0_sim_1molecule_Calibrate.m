
% Script to simulate snapshots of molecules using simSPT trajectory and
% save noisy snapshots to calibrate conversion factor between PT area and diffusivity.

%% DESCRIPTION
% 1) Simulate snapshot of molecules using simSPT trajectory (simulation
% time gap 1ms). One snapshot is generated from one trajectory with
% desired track length from simSPT. The pixelized image pair
% have the exactly 0.11um/pix size.
% 2) A Gaussian noise and Possion modeled shot-noise is added into simulated snapshot.

% Citation of Hansen, Anders S., Iryna Pustova, Claudia Cattoglio, Robert Tjian, and Xavier Darzacq. “CTCF and Cohesin Regulate Chromatin Loop Stability with Distinct Dynamics.” ELife 6 (May 3, 2017): e25776. https://doi.org/10.7554/eLife.25776. 
% simSPT: "https://gitlab.com/tjiandarzacq-lab/simSPT"

%% Zuhui Wang
%% 2023/06/09
%% INTENTIONALLY KEEP BLANK
%% INTENTIONALLY KEEP BLANK
%% INTENTIONALLY KEEP BLANK

%% Obtain real immobile PSF from experimental data
close all; clear; rng(37, 'twister'); % for reproducibility
[folderPath,~,~] = fileparts(mfilename('fullpath'));
addpath(genpath(folderPath)) % add folder of step1 to the search path

% change your desired imaging parameters
target_pixelSize = 0.11; % desired simulated camera pixel size 0.11 or 0.16, um
background_mean = [1200 1200]; % mean value of background, unit ADU, can obtain from ImageJ by measuring the mean value of background region
background_std = [300 300];% std value of background, unit ADU, can obtain from ImageJ by measuring the std value of background region
good_snr_range = [35 35]; % SNR range for molecule detections, unit dB, for calibration, set to high SNR range 35
background_offset = 500; % background offset of camera when no incident photon on camera
num_single_images_to_generate = 100; % number of single images to generate under each diffusivity, 100 is enough
PhotonConversionRate = 0.9398923112; % (0.5e/ADU, EMGain300, QE at 660nm:0.9398923112) (0.5e/ADU, EMGain300, QE at 582nm:0.9696909185)

% Please follow the extract_psf.ipynb to obtain your PSF sigma value
beads_PixelSize=0.11; % um per pixel, pixel size of beads image
mean_r0 = 1.4141*beads_PixelSize; % replace 1.4141 with your microscope's sigma value

% Reference:
% simSPT: "https://gitlab.com/tjiandarzacq-lab/simSPT" from Hansen, Anders S., Iryna Pustova, Claudia Cattoglio, Robert Tjian, and Xavier Darzacq. “CTCF and Cohesin Regulate Chromatin Loop Stability with Distinct Dynamics.” ELife 6 (May 3, 2017): e25776. https://doi.org/10.7554/eLife.25776. 
% msdanalyzer: Jean-Yves Tinevez (2025). Mean square displacement analysis of particles trajectories (https://github.com/tinevez/msdanalyzer), GitHub. Retrieved January 11, 2025.


%% Simulate snapshots from simSPT
% import simSPT with time delay 1ms
simSPT_path = '/path/to/simSPT_generated_tracks';
simSPT_file = [...
    % >>>>>>>> No localization error >>>>>>>
    "20231202_logDm3p0_pure_10kTraj_dt1ms.csv",...
    "20231202_logDm2p8_pure_10kTraj_dt1ms.csv",...
    "20231202_logDm2p6_pure_10kTraj_dt1ms.csv",...
    "20231202_logDm2p4_pure_10kTraj_dt1ms.csv",...
    "20231202_logDm2p2_pure_10kTraj_dt1ms.csv",...
    "20231202_logDm2p0_pure_10kTraj_dt1ms.csv",...
    "20231202_logDm1p8_pure_10kTraj_dt1ms.csv",...
    "20231202_logDm1p6_pure_10kTraj_dt1ms.csv",...
    "20231202_logDm1p4_pure_10kTraj_dt1ms.csv",...
    "20231202_logDm1p2_pure_10kTraj_dt1ms.csv",...
    "20231202_logDm1p0_pure_10kTraj_dt1ms.csv",...
    "20231202_logDm0p8_pure_10kTraj_dt1ms.csv",...
    "20231202_logDm0p6_pure_10kTraj_dt1ms.csv",...
    "20231202_logDm0p4_pure_10kTraj_dt1ms.csv",...
    "20231202_logDm0p2_pure_10kTraj_dt1ms.csv",...
    "20231202_logD0p0_pure_10kTraj_dt1ms.csv",...
    "20231202_logD0p2_pure_10kTraj_dt1ms.csv",...
    "20231202_logD0p4_pure_10kTraj_dt1ms.csv",...
    "20231202_logD0p6_pure_10kTraj_dt1ms.csv",...
    "20231202_logD0p8_pure_30kTraj_dt1ms.csv",...    
    "20231202_logD1p0_pure_30kTraj_dt1ms.csv",...         
    % <<<<<<<<< No localization error <<<<<<<<<<
    ];

simSPT_D = [...
    "logDm3p0",...
    "logDm2p8",...
    "logDm2p6",...
    "logDm2p4",...
    "logDm2p2",...
    "logDm2p0",...
    "logDm1p8",...
    "logDm1p6",...
    "logDm1p4",...
    "logDm1p2",...
    "logDm1p0",...
    "logDm0p8",...
    "logDm0p6",...
    "logDm0p4",...
    "logDm0p2",...
    "logD0p0",...
    "logD0p2",...
    "logD0p4",...
    "logD0p6",...
    "logD0p8",...    
    "logD1p0",...    
    ];

DiffCoeff = [0.00100    0.00158	0.00251	0.00398	0.00631	0.01000	0.01585	0.02512	0.03981	0.06310	0.10000	0.15849	0.25119	0.39811	0.63096	1.00000	1.58489	2.51189	3.98107	6.30957	7.94328	10.00000];

% desired exposure time to generate PSF
target_exposure = 30; % unit ms

windowSize_array = [ones(1,10)*50,ones(1,6)*10,ones(1,6)*5];
max_track_array = [ones(1,22)*100];

% parent directory to save image pairs
img_pair_saveDir = '/path/to/save/simulated_dataset_forCalibrate_DiffCoeff'; % please change to your desired path
mkdir(fullfile(img_pair_saveDir));
mkdir(fullfile(img_pair_saveDir,'imgs'));
mkdir(fullfile(img_pair_saveDir,'annotations'));


%% ================ Start snapshots simulation ================ %%
% ----------------- Prepare track list ----------------------- %
all_tracks = cell(1,length(simSPT_file));
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
        windowSize = windowSize_array(iSPT);
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
                        if total_subtrack_idx > max_track_array(iSPT)
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

        fprintf('Total individual trajectories prepared: %d\n', length(long_track_id));
        
        now_track_cell = cell(1, length(long_track_id));
        for k = 1:length(long_track_id)
            now_track_cell{k} = long_t(long_t.trajecotry_newIdx == long_track_id(k),:);
        end
        
        if length(now_track_cell) < max_track_array(iSPT)
            fprintf('Current trajectory# in this condition %s is lower than max_imgpair!',simSPT_file{iSPT});
        end

        all_tracks{iSPT} = now_track_cell;
               
    end % End loop over target_exposure    
end % End loop over simSPT_file

%% ----------------- Generate multi-molecule images ---------------------- %
fprintf('Will generate %d multi-molecule images.\n', num_single_images_to_generate);

% Loop to generate multiple multi-molecule images

num_simSPT_file = length(simSPT_file);
num_tracks = cellfun(@length,all_tracks);

% NumWorkers = 8;
% parpool('local', NumWorkers);
for iSPT = 1:length(simSPT_file)
    tic;
    fileID = fopen(fullfile(img_pair_saveDir,'annotations',sprintf('sample_D%s.info.txt',simSPT_D(iSPT))),'w');
    fprintf(fileID,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', ...
        'frame','loc_id','simSPT_D','D_value','D_fromMSD','logD_fromMSD','b_fromMSD','r2_fromMSD','exposure_time','molecule_class','x_pix_pos','y_pix_pos','SNR','background','noise');
    progressbarText(0)
    for img_count = 0:num_single_images_to_generate-1
        progressbarText(img_count/num_single_images_to_generate);
        mol_count = 1;   
        iTrack = randi([1, num_tracks(iSPT)]);                    
        now_track = all_tracks{iSPT}{iTrack};            
        
        % simulate a snapshots by merge individual fitted psf from the same
        % trajectory

        % Calculate the range of x and y to span 3.52um (3.52/0.11=32
        % pixels) so that later can bin into 32x32 pixel with pixel
        % size 0.11um/px

        % Calculate the range of x and y to span 5.12um (5.12/0.16=32
        % pixels) so that later can bin into 32x32 pixel with pixel
        % size 0.16um/px

        range = 32*target_pixelSize/2;
        interval = 0.01; % spacing between x and y values
        x_min = mean(now_track.x) - range;
        x_max = mean(now_track.x) + range;
        y_min = mean(now_track.y) - range;
        y_max = mean(now_track.y) + range;
        
        % Generate x and y vectors that span 2xrange with 0.01 um gap size
        x = x_min:interval:x_max;
        y = y_min:interval:y_max;

        [Xgrid,Ygrid] = meshgrid(x,y);
    
        f = @(I,r0,m,x0,y0,x,y) I * exp(-(1/(2*r0^2))*((y - y0).^2 + (x - x0).^2)) + m;
    
        merge_z = zeros(length(Ygrid),length(Xgrid));            

        for i = 1:height(now_track)                                
            rand_I = 1;
            rand_r0 = mean_r0;
            rand_m = 0; % no background consider the ideal situation        
            z = feval(f,rand_I,rand_r0,rand_m,now_track.x(i),now_track.y(i),Xgrid,Ygrid);                
            merge_z = merge_z + z;        
        end

        merge_z = merge_z/sum(merge_z(:));  % Normalize the matrix by dividing by the total volume
        
        % Grid the snapshotsry image to match the real camera pixel
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
        trackImg_pixelSize = target_pixelSize/10; % unit um; pixel size of finer trajectory image
        trackImg_range = 320*trackImg_pixelSize/2;
        x_min = mean(now_track.x) - trackImg_range;
        x_max = mean(now_track.x) + trackImg_range;
        y_min = mean(now_track.y) - trackImg_range;
        y_max = mean(now_track.y) + trackImg_range;
        
        pix_track_x = round((now_track.x - x_min) / trackImg_pixelSize + 0.5);
        pix_track_y = round((now_track.y - y_min) / trackImg_pixelSize + 0.5);

        OutOfRangeError = true;
        while OutOfRangeError
            % Check if the track is within the 320x320 image range                
            OutOfRangeError = false;
            for i = 2:length(pix_track_x)
                [x_line, y_line] = bresenham(pix_track_x(i-1), pix_track_y(i-1), pix_track_x(i), pix_track_y(i));
                if any(x_line < 1) || any(x_line > 320) ||  any(y_line < 1) || any(y_line > 320)
                    OutOfRangeError = true;
                    break; % Current track cannot fit 320x320 image
                end
            end
            
            
            if min(pix_track_x) < 0 || min(pix_track_y) < 0 || max(pix_track_x) > 320 || max(pix_track_y) > 320 || OutOfRangeError

                continue % Current track cannot fit 320x320 image

            else
                            
                %=========    create ground truth of motion blur with good SNR ===============%                                                               
                % Set the range of the uniform distribution of SNR
                SNR_range = good_snr_range; %[17 30]; 
                % Generate a single random number from the uniform distribution
                SNR = SNR_range(1) + (SNR_range(2)-SNR_range(1))*rand(1,1);
                bk_m = background_mean(1) + (background_mean(2)-background_mean(1))*rand(1,1);% background_mean;%1761;% 1152;
                bk_sig = background_std(1) + (background_std(2)-background_std(1))*rand(1,1);%background_std;%409;% 283;                
                rescale_pix_merge_z = rescale(pix_merge_z);                

                % Noise model : Poisson shot-noise + Gaussian white noise
                noisy_pix_merge_z = uint16(poissrnd(rescale_pix_merge_z*10^((SNR+20*log10(bk_sig))/20)) + bk_m + bk_sig*randn(32,32));
                noisy_pix_merge_z(noisy_pix_merge_z<background_offset) = background_offset;

                % Initialize the msdanalyzer object
                temp_simSPT_track_msd = msdanalyzer(2,'um','seconds');                      
                temp_simSPT_track_msd = temp_simSPT_track_msd.addAll({[now_track.t,now_track.x,now_track.y]});                
                % calculate the MSD of all imported tracks
                temp_simSPT_track_msd = temp_simSPT_track_msd.computeMSD([]);
                % calculate D from MSD-dT, use 0.25 clip factor, use origin as fitting point including weight
                temp_simSPT_track_msd = temp_simSPT_track_msd.fitMSD();
                                
                D_fromMSD = temp_simSPT_track_msd.lfit.a/4; % D_fromMSD
                b_fromMSD = temp_simSPT_track_msd.lfit.b; % b_fromMSD
                r2_fromMSD = temp_simSPT_track_msd.lfit.r2fit; % r2_fromMSD
                logD_fromMSD = log10(D_fromMSD);

                placed_center_x_px = 16.5-1; % start from 0
                placed_center_y_px = 16.5-1; % start from 0

                fprintf(fileID,'%8d\t%8d\t%8s\t%8.3f\t%8.3f\t%8.6f\t%8.6f\t%8.6f\t%3d\t%s\t%8.3f\t%8.3f\t%8.3f\t%8.3f\t%8.3f\n', ...                    
                    img_count,mol_count,simSPT_D(iSPT),DiffCoeff(iSPT),D_fromMSD,logD_fromMSD,b_fromMSD,r2_fromMSD,target_exposure,'molecule',placed_center_x_px,placed_center_y_px,SNR,bk_m,bk_sig);                
                                    
                if img_count == 0                
                    % Save the multi-molecule image
                    outfile_img = fullfile(img_pair_saveDir,'imgs', ...
                        sprintf('sample_D%s.img.tif', simSPT_D(iSPT)));            
                    imwrite(uint16(noisy_pix_merge_z), outfile_img);                            
                else
                    imwrite(uint16(noisy_pix_merge_z), outfile_img, 'WriteMode', 'append');                         
                end                                                                     
            end
        end                         
    end
    fclose(fileID); 
    toc;  
end


% copy this script to the data folder as a log
sourceFile = mfilename('fullpath'); % Define the source and target file paths
copyfile([sourceFile '.m'], img_pair_saveDir);

disp("Simulated ground truth export finished!");
