% Script to simulate motion-blur of molecules using simSPT trajectory and
% save noisy motion-blur and binary mask molecule image pairs

%% DESCRIPTION
% 1) Simulate motion-blur of molecules using simSPT trajectory (simulation
% time gap 1ms). One motion-blur is generated from one trajectory with
% desired track length from simSPT. Then, save the motion-blur and sharp
% (immobile) molecule image pair into one image and will be used as the
% ground truth of training Motion-ETR network. The pixelized image pair
% have the exactly 0.11um/pix size.
% 2) A Gaussian noise is added into simulated motion-blur.
% 3) Simulated motion-blur is further do ellipticity index fitting.

% Careful with the definition of PSF parameters (I, r0, m) For r0, I use impars.psfStd from MTT.
% Dependency: bfmatlab >6.6.0 

% Citation of Hansen, Anders S., Iryna Pustova, Claudia Cattoglio, Robert Tjian, and Xavier Darzacq. “CTCF and Cohesin Regulate Chromatin Loop Stability with Distinct Dynamics.” ELife 6 (May 3, 2017): e25776. https://doi.org/10.7554/eLife.25776. 
% simSPT: "https://gitlab.com/tjiandarzacq-lab/simSPT"

%% Zuhui Wang
%% 2023/06/09
%% INTENTIONALLY KEEP BLANK
%% INTENTIONALLY KEEP BLANK
%% INTENTIONALLY KEEP BLANK



%% Obtain real immobile PSF from experimental data
close all; clear;clc;
addpath(genpath('/home/zuhui2/Documents/METHOD/simulation/simPSF'),'-end')
addpath(genpath('/home/zuhui2/Documents/MATLAB/msdanalyzer'),'-end')

%%% IMPORTANT CURRENTLY ALL PARAMS SET TO SIMULATE 0.11UM PIXEL SIZE
% change your desired imaging parameters
target_pixelSize = 0.11;
EmissionWavelength = 664;%664; % 580 for mEos3.2 & mEosEM
save_MSD = false; % save MSD of each trajectories

% import raw image and psf data
img_path = '/dataB/zuhui2/Simulation/MotionBlurDataset/';
img_name = '20221114_Cell01_FOXA2-Halo_10uMPA646_10p5ms_2kframe_01.nd2';
imgs_3d_matrix = MemoryEfficientND2reader(fullfile(img_path,img_name));

psf_path = '/dataE/WZH-DataCenter/PROCESS-SPT/2022/20221109-10_PA646_U2OS_Xlone_FOXA2-Halo_testPSFExposurePixelSize/SNR_15/20221114_Cell01_FOXA2-Halo_10uMPA646_10p5ms_2kframe_01/';
psf_name = 'v3_5_singleCell_2022-11-16_23_23_33_PSFSNR15_COVp25_PSF.mat';
psffit_name = 'v3_5_singleCell_2022-11-16_23_23_33_PSFSNR15_COVp25_PSF_FitResult.csv';


load(fullfile(psf_path,psf_name),'cell_PSF')
fittable = readtable(fullfile(psf_path,psffit_name));
impars.PixelSize=target_pixelSize; % um per pixel
impars.psf_scale=1.35; % PSF scaling
impars.wvlnth= EmissionWavelength/1000; %emission wavelength in um
impars.NA=1.49; % NA of detection objective
impars.psfStd= impars.psf_scale*0.55*(impars.wvlnth)/impars.NA/1.17/impars.PixelSize/2; % PSF standard deviation in pixels

cell_PSF_xyI = vertcat(cell_PSF.xyI{:});
cell_PSF_BoundxyI = vertcat(cell_PSF.BoundxyI{:});

% remove low quality PSF
PSFSz = cellfun(@height,cell_PSF_xyI);
fittable.PSF_Sz = PSFSz;
fittable.SNRvSz = fittable.SNR./fittable.PSF_Sz;
cell_PSF_xyI = cell_PSF_xyI(fittable.SNRvSz > 0.25);
cell_PSF_BoundxyI = cell_PSF_BoundxyI(fittable.SNRvSz > 0.25);
fittable = fittable(fittable.SNRvSz > 0.25,:);

% obtain immobile psf
immobile_psfIdx = find(fittable.EllipticityIdx<1);
immo_cell_PSF_xyI = cell_PSF_xyI(immobile_psfIdx);
immo_cell_PSF_BoundxyI = cell_PSF_BoundxyI(immobile_psfIdx);

% % check some psf
% for i =  randi(length(immobile_psfIdx),1,10)
%     PSF_3DIMG(fittable.Frame(immobile_psfIdx(i)),fittable.PSF_idx(immobile_psfIdx(i)),cell_PSF)
% end

% fit immobile psf
[Xgrid,Ygrid,Zgrid] = cellfun(@ARRAY2GRID,immo_cell_PSF_xyI,immo_cell_PSF_BoundxyI,'UniformOutput',false);
[xData, yData, zData] = cellfun(@prepareSurfaceData,Xgrid,Ygrid,Zgrid,'UniformOutput',false);

cell_fitresult = cellfun(@cell_fit_gaussian2,xData,yData,zData,'UniformOutput',false);
cell_fitresult = vertcat(cell_fitresult{:});

% simulate single psf using the fit psf function with added image noise
mean_I = mean(cell_fitresult(:,1));
sigma_I = std(cell_fitresult(:,1));
% mean_r0 = mean(cell_fitresult(:,2)); 
% sigma_r0 = std(cell_fitresult(:,2));
mean_r0 = impars.psfStd*impars.PixelSize; % unit in um, use MTT calculated PSF std
sigma_r0 = 0;% unit in um, use MTT calculated PSF std, no variation
% image offset estimation from fitting
mean_m = mean(cell_fitresult(:,3));
sigma_m = std(cell_fitresult(:,3)); % image noise estimation from raw image


% % (debug only) visualize the simulated PSF
% f = @(I,r0,m,x0,y0,x,y)(I/r0) * exp(-(1/(2*r0^2))*((y - y0).^2 + (x - x0).^2)) + m;
% x = 1:15;
% y = 1:15;
% [Xgrid,Ygrid] = meshgrid(x,y);
% noise_z = normrnd(0,sigma_m,15,15);
% rand_I = normrnd(mean_I,sigma_I);
% rand_r0 = normrnd(mean_r0,sigma_r0);
% rand_m = normrnd(mean_m,sigma_m);
% offset_addnoise = rand_m + noise_z;
% while (rand_I < 0 || rand_r0  < 0 || any(offset_addnoise(:) < 0))
%     noise_z = normrnd(0,sigma_noise,15,15);
%     rand_I = normrnd(mean_I/10,sigma_I/10);
%     rand_r0 = normrnd(mean_r0,sigma_r0/10);
%     rand_m = normrnd(mean_m,sigma_m);
%     offset_addnoise = rand_m + noise_z;
% end
% z = feval(f,rand_I,rand_r0,rand_m,mean(Xgrid(:)),mean(Ygrid(:)),Xgrid,Ygrid);
% figure;surf(z+noise_z)


%% Simulate motion-blur from simSPT
% import simSPT with time delay 1ms
% simSPT code: ./simSPT -D1=0.03 -D2=2.0 -p1=0.5 -p2=0.5 -sigma=0.031 -dt=0.001 -n_traj=100000 -file=../simPSF/20230309_D2_p5_dt1ms.csv -seed=0
simSPT_path = '/home/zuhui2/Documents/METHOD/simulation/simPSF/Results/20231204_LocError60nm';
simSPT_file = [...
    % "20230314/20230314_D01_pure_10kTraj_dt1ms.csv",...
    % "20230314/20230314_D02_pure_10kTraj_dt1ms.csv",...
    % "20230314/20230314_D03_pure_10kTraj_dt1ms.csv",...
    % "20230314/20230314_D04_pure_10kTraj_dt1ms.csv",...
    % "20230314/20230314_D05_pure_10kTraj_dt1ms.csv",...
    % "20230314/20230314_D06_pure_10kTraj_dt1ms.csv",...
    % "20230314/20230314_D07_pure_10kTraj_dt1ms.csv",...
    % "20230314/20230314_D08_pure_10kTraj_dt1ms.csv",...
    % "20230314/20230314_D09_pure_10kTraj_dt1ms.csv",...
    % "20230314/20230314_D10_pure_10kTraj_dt1ms.csv",...
    % "20230606/20230606_D0p01_pure_10kTraj_dt1ms.csv",...
    % "20230606/20230606_D0p1_pure_10kTraj_dt1ms.csv",...
    % "20230606/20230606_D01_pure_10kTraj_dt1ms.csv",...
    % "20230606/20230606_D10_pure_10kTraj_dt1ms.csv",...
    % % >>>>>>>> With 31 nm localization error >>>>>>>
    % "20230612_D0p01_pure_10kTraj_dt1ms.csv",...
    % "20230612_D0p1_pure_10kTraj_dt1ms.csv",...
    % "20230612_D01_pure_10kTraj_dt1ms.csv",...
    % "20230612_D02_pure_10kTraj_dt1ms.csv",...
    % "20230612_D03_pure_10kTraj_dt1ms.csv",...
    % "20230612_D04_pure_10kTraj_dt1ms.csv",...
    % "20230612_D05_pure_10kTraj_dt1ms.csv",...
    % "20230612_D06_pure_30kTraj_dt1ms.csv",...
    % "20230612_D07_pure_30kTraj_dt1ms.csv",...
    % "20230612_D08_pure_30kTraj_dt1ms.csv",...
    % "20230612_D09_pure_30kTraj_dt1ms.csv",...
    % "20230613_D10_pure_35kTraj_dt1ms.csv",...
    % "20230612_D15_pure_10kTraj_dt1ms.csv",...
    % "20230612_D20_pure_10kTraj_dt1ms.csv",...
    % "20230612_D25_pure_10kTraj_dt1ms.csv",...
    % "20230612_D30_pure_10kTraj_dt1ms.csv",...
    % "20230612_D35_pure_10kTraj_dt1ms.csv",...
    % % <<<<<<<<< With 31 nm localization error <<<<<<<<<<
    % % >>>>>>>> No localization error >>>>>>>
    % "20231202_D0p01_pure_10kTraj_dt1ms.csv",...
    % "20231202_D0p1_pure_10kTraj_dt1ms.csv",...
    % "20231202_D01_pure_10kTraj_dt1ms.csv",...
    % "20231202_D02_pure_10kTraj_dt1ms.csv",...
    % "20231202_D03_pure_10kTraj_dt1ms.csv",...
    % "20231202_D04_pure_10kTraj_dt1ms.csv",...
    % "20231202_D05_pure_10kTraj_dt1ms.csv",...
    % "20231202_D06_pure_30kTraj_dt1ms.csv",...
    % "20231202_D07_pure_30kTraj_dt1ms.csv",...
    % "20231202_D08_pure_30kTraj_dt1ms.csv",...
    % "20231202_D09_pure_30kTraj_dt1ms.csv",...
    % "20231202_D10_pure_10kTraj_dt1ms.csv",...
    % "20231202_D15_pure_10kTraj_dt1ms.csv",...
    % "20231202_D20_pure_10kTraj_dt1ms.csv",...
    % "20231202_D25_pure_10kTraj_dt1ms.csv",...
    % "20231202_D30_pure_10kTraj_dt1ms.csv",...
    % "20231202_D35_pure_10kTraj_dt1ms.csv",...
    % % <<<<<<<<< No localization error <<<<<<<<<<
    % >>>>>>>> With 60 nm localization error >>>>>>>
    "20231204_D0p01_pure_10kTraj_dt1ms.csv",...
    "20231204_D0p1_pure_10kTraj_dt1ms.csv",...
    "20231204_D01_pure_10kTraj_dt1ms.csv",...
    "20231204_D02_pure_10kTraj_dt1ms.csv",...
    "20231204_D03_pure_10kTraj_dt1ms.csv",...
    "20231204_D04_pure_10kTraj_dt1ms.csv",...
    "20231204_D05_pure_10kTraj_dt1ms.csv",...
    "20231204_D06_pure_30kTraj_dt1ms.csv",...
    "20231204_D07_pure_30kTraj_dt1ms.csv",...
    "20231204_D08_pure_30kTraj_dt1ms.csv",...
    "20231204_D09_pure_30kTraj_dt1ms.csv",...
    "20231204_D10_pure_10kTraj_dt1ms.csv",...
    "20231204_D15_pure_10kTraj_dt1ms.csv",...
    "20231204_D20_pure_10kTraj_dt1ms.csv",...
    "20231204_D25_pure_10kTraj_dt1ms.csv",...
    "20231204_D30_pure_10kTraj_dt1ms.csv",...
    "20231204_D35_pure_10kTraj_dt1ms.csv",...
    % <<<<<<<<< With 60 nm localization error <<<<<<<<<<
    ];

% corresponding DiffusionCoefficient to the imported simSPT_file
simSPT_D_num = [0.01 0.1 1:10 15 20 25 30 35]; 
simSPT_D = ["p01", "p1", "01","02","03","04","05","06","07","08","09","10","15","20","25","30","35"]; 
% Starting index of simSPT trajectories
StartIdx = 1;
% maximal number of generated image pair at given D and exposure time
max_imgpair = 5000; 
% desired exposure time to generate PSF
target_exposure = 1;%30; %[1:10 20:10:60]; 

% parent directory to save image pairs
img_pair_saveDir = '/dataE/WZH-DataCenter/PROCESS-SPT/2023/simPSF_results/UNetDetectionDataset/20231204_NoNoiseLocError60_1msTimeGap_110nmPixSize/D_upto35/';
mkdir(img_pair_saveDir)
% mkdir(fullfile(img_pair_saveDir,'imgs'));
% mkdir(fullfile(img_pair_saveDir,'masks'));

% Table to save area of smooth or pixelized area of simulated motion-blur
SummaryTable_savePath = img_pair_saveDir;
if save_MSD
    SummaryTable = table('Size',[0,8],'VariableTypes',{'cell','double','double','double','double','double','double','double'},...
        'VariableNames',{'ImageName','new_track_id','ExposureTime','DiffCoeff','MSD','SmoothArea','PixelizedArea','JumpStartToEnd'}); 
    writetable(SummaryTable,fullfile(SummaryTable_savePath,sprintf('NoNoise_SmoothArea_PixelArea.csv')),'WriteRowNames',true);
else
    SummaryTable = table('Size',[0,7],'VariableTypes',{'cell','double','double','double','double','double','double'},...
        'VariableNames',{'ImageName','new_track_id','ExposureTime','DiffCoeff','SmoothArea','PixelizedArea','JumpStartToEnd'}); 
    writetable(SummaryTable,fullfile(SummaryTable_savePath,sprintf('NoNoise_SmoothArea_PixelArea.csv')),'WriteRowNames',true);
end

% initiate the log file to save image pair info
fileID = fopen(fullfile(img_pair_saveDir,'Summary.txt'),'w');
fprintf(fileID,'%s\t%g\n','target_pixelSize',target_pixelSize);
fprintf(fileID,'%s\t%g\n','EmissionWavelength',EmissionWavelength);
% fprintf(fileID,'%s\t%d\n','background_mean',background_mean);
% fprintf(fileID,'%s\t%d\n','background_std',background_std);
% fprintf(fileID,'%s\t%d\t%d\n','good_snr_range',good_snr_range(1),good_snr_range(2));
% fprintf(fileID,'%s\t%d\t%d\n','bad_snr_range',bad_snr_range(1),bad_snr_range(2));
fprintf(fileID,'%s\t%s\t%s\n','simSPT_D','exposureTime','imagePair#');

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
        % Note: here we want to simulate frame time gap 1ms, not exposure 1ms, so I plus 1 to get 2 localizations per track
        windowSize = 1;
        long_t = table();
        stopAllLoops = false; % Initialize flag variable
        total_subtrack_idx = 1;
        for i = 1:length(track_id)
            trajIndex = track_id(i);
            trajData = t_raw(t_raw.trajectory == trajIndex, :);
            dataHeight = height(trajData);

            if dataHeight >= target_exposure(idT)+1
                for j = 1:windowSize:dataHeight
                    if j+target_exposure(idT) <= dataHeight % not exceed the data height
                        subTrajData = trajData(j:j+target_exposure(idT), :);
                        subTrajData.trajecotry_newIdx = ones(target_exposure(idT)+1,1)*total_subtrack_idx;
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
        progressbarText(0);
        for iter = 1:length(long_track_id)
            progressbarText(iter/length(long_track_id));
            now_track_cell{iter} = long_t(long_t.trajecotry_newIdx == long_track_id(iter),:);
        end

        if length(now_track_cell) < max_imgpair
            fprintf('Current trajectory# in this file is lower than max_imgpair!');
        end
        
        % initiate variable
        pix_Xgrid_cell = {};
        pix_Ygrid_cell = {};
        pix_mergeZ_cell = {};
        pix_masked_mergeZ_cell = {};
        pix_filename = {};
        mask_pix_merge_z_area = nan(length(now_track_cell),1);
        mask_merge_z_area = nan(length(now_track_cell),1);
        JumpStartToEnd = nan(length(now_track_cell),1);

        % parfor_progress(length(now_track_cell));
        parfor iter = 1:length(now_track_cell)
            %progressbarText(iTrack/length(long_track_id));
            now_track = now_track_cell{iter};
           
            % simulate a motion-blur by merge individual fitted psf from the same
            % trajectory

            % Calculate the range of x and y to span 3.52um (3.52/0.11=32
            % pixels) so that later can bin into 32x32 pixel with pixel
            % size 0.11um/px
            % if target_pixelSize == 0.11
            %     range = 1.76; % range of x and y from center position; 1.76 for 0.11um pixel size; 2.56 for 0.16um pixel size
            % elseif target_pixelSize == 0.16
            %     range = 2.56;
            % else
            %     f = errordlg('target_pixelSize not right','Error');
            % end
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
        
            f = @(I,r0,m,x0,y0,x,y)(I/r0) * exp(-(1/(2*r0^2))*((y - y0).^2 + (x - x0).^2)) + m;
        
            merge_z = zeros(length(Ygrid),length(Xgrid));
            tempSharp_z = zeros(length(Ygrid),length(Xgrid),height(now_track)); % store instant frame intensity

            for i = 1:height(now_track)
                
                % noise_z = zeros(length(Ygrid),length(Xgrid)); % no noise consider the ideal situation
        
                rand_I = normrnd(mean_I/10,sigma_I/100);
                % rand_r0 = normrnd(mean_r0,sigma_r0);
                rand_r0 = mean_r0;
                rand_m = 0; % no background consider the ideal situation
        
                while (rand_I < 0 || rand_r0  < 0 || rand_m < 0 )
                    rand_I = normrnd(mean_I/10,sigma_I/100);
                    % rand_r0 = normrnd(mean_r0,sigma_r0);
                    rand_r0 = mean_r0;
                    rand_m = 0; 
                end
        
                z = feval(f,rand_I,rand_r0,rand_m,now_track.x(i),now_track.y(i),Xgrid,Ygrid);
                % tempSharp_z(:,:,i) = feval(f,rand_I,rand_r0,rand_m,mean(now_track.x),mean(now_track.y),Xgrid,Ygrid); % sharp image
                merge_z = merge_z + z;
        
            end

            merge_z = merge_z/sum(merge_z(:));  % Normalize the matrix by dividing by the total volume
            
            % % create sharp image with averaged pixel intensity (matrix z)
            % mean_sharp_z = mean(tempSharp_z,3);

            % Grid the motion-blurry image to match the real camera pixel
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

            % ----------- GENERATE MASK ------------%
            % generate the pixelized image pair, rescale the merged image to the range [0, 1]
            % Create binary mask of motion blur, maxentropie only takes uint8 image
            % [~, mask_pix_merge_z]=maxentropie(im2uint8(rescale(pix_merge_z)));
            % mask_pix_merge_z_area(iter) = (sum(mask_pix_merge_z,'all')/255)*target_pixelSize*target_pixelSize;
            % mask_pix_merge_z = im2uint8(rescale(mask_pix_merge_z)); % convert from double to uint8
            
            % mask_merge_z = imbinarize(rescale(merge_z),0.14);
            % mask_merge_z_area(iter) = sum(mask_merge_z,'all')*interval*interval;
            % mask_merge_z = im2uint8(rescale(mask_merge_z)); % convert from double to uint8

            % maxentropie equally cover >95% of volumn by testing on
            % immobile simulated molecule

            % Sort the normalized volume values in descending order
            sorted_values = sort(merge_z(:), 'descend');

            % Find the index where the cumulative sum exceeds 99% of the total volume
            cumulative_sum = cumsum(sorted_values);
            threshold_index = find(cumulative_sum >= 0.95, 1);

            % Create a binary mask using the threshold index
            mask_merge_z = merge_z >= sorted_values(threshold_index);
            mask_merge_z_area(iter) = sum(mask_merge_z,'all')*interval*interval;

            % Sort the normalized volume values in descending order
            sorted_values = sort(pix_merge_z(:), 'descend');

            % Find the index where the cumulative sum exceeds 99% of the total volume
            cumulative_sum = cumsum(sorted_values);
            threshold_index = find(cumulative_sum >= 0.95, 1);

            % Create a binary mask using the threshold index
            mask_pix_merge_z = pix_merge_z >= sorted_values(threshold_index);
            mask_pix_merge_z_area(iter) = sum(mask_pix_merge_z,'all')*target_pixelSize*target_pixelSize;

            
            % save pixel image into cell for PSF fitting
            % calculate grid center
            % pix_center_x = pix_Xedge(1:end-1)+(pix_Xedge(2)-pix_Xedge(1))/2;
            % pix_center_y = pix_Yedge(1:end-1)+(pix_Yedge(2)-pix_Yedge(1))/2;         
            % [pix_Xgrid, pix_Ygrid] = meshgrid(pix_center_x,pix_center_y);
            % pix_Xgrid_cell{iter} = pix_Xgrid;
            % pix_Ygrid_cell{iter} = pix_Ygrid;
            % pix_mergeZ_cell{iter} = pix_merge_z;

            % save pixel image into cell for PSF fitting
            smooth_Xgrid_cell{iter} = Xgrid;
            smooth_Ygrid_cell{iter} = Ygrid;
            smooth_mergeZ_cell{iter} = merge_z;

            JumpStartToEnd(iter) = sqrt((now_track.x(end) - now_track.x(1))^2 + (now_track.y(end) - now_track.y(1))^2);

            % % (Debug only) Visualize the smooth and pixelized motion-blur and its trajectory
            % figure;imshow(rescale(merge_z))
            % figure;imshow(mask_merge_z)
            % figure;imshow(mask_pix_merge_z,'InitialMagnification',1000);
            % figure;imshow(rescale(pix_merge_z),'InitialMagnification',1000)

            % figure;plot(now_track.x,now_track.y,'LineWidth',1)
            % xlim([x_min,x_max])
            % ylim([y_min,y_max])
            % axis ij equal
            
            % imwrite(pix_merge_z,fullfile(img_pair_saveDir,'imgs',sprintf('%s_TrackID%012d_D%s_dT%02dms.tif','Pixelized',now_track.trajecotry_newIdx(1),simSPT_D(iSPT),target_exposure(idT))));
            % imwrite(mask_pix_merge_z,fullfile(img_pair_saveDir,'masks',sprintf('%s_TrackID%012d_D%s_dT%02dms_mask.gif','Pixelized',now_track.trajecotry_newIdx(1),simSPT_D(iSPT),target_exposure(idT))));

            % File Name
            pix_filename{iter} = sprintf('%s_TrackID%012d_D%s_dT%02dms','Pixelized',now_track.trajecotry_newIdx(1),simSPT_D(iSPT),target_exposure(idT));
            
        end
        

        % %% ------------- PSF fitting ------------------%%
        % fprintf('PSF fitting ...\n');    
        % % [xData, yData,zData] = cellfun(@prepareSurfaceData,smooth_Xgrid_cell,smooth_Ygrid_cell, smooth_mergeZ_cell, 'UniformOutput', false);
        
        % % Version 1 CPU Parallel Computing-fast to do PSF fitting
        % % ATTENTION: CHANGE make_BoundxyI INDEX WHEN PIXEL IMAGE DIMENSION
        % % CHANGED
        % % sim_PSF_xyI = cellfun(@horzcat,xData, yData, zData, 'UniformOutput', false);
        % % sim_PSF_BoundxyI = cellfun(@make_BoundxyI,smooth_Xgrid_cell,smooth_Ygrid_cell,smooth_mergeZ_cell,'UniformOutput', false);
        
        % % fitresult = PARA_GAUSS2DFIT(cell_PSF_xyI,cell_PSF_BoundxyI);
        % % [fitresult, ~, ~]= PARA_GAUSS2DFIT(sim_PSF_xyI',sim_PSF_BoundxyI');
        % % [fitresult, ~, ~]= PARA_SIM_GAUSS2DFIT(sim_PSF_xyI');
        % [fitresult, ~, ~]= PARA_SIM_GAUSS2DFIT(smooth_Xgrid_cell,smooth_Ygrid_cell, smooth_mergeZ_cell);
        
        % % Append ellipticity
        % fitresult(:,8) = max(fitresult(:,3:4),[],2)./min(fitresult(:,3:4),[],2);
        % % Append ellipticity index
        % fitresult(:,9) = 100*log10(fitresult(:,8));

        %% ------------------- calculate MSD of simSPT trajectories ------------------- %%
        if save_MSD
            % Initialize the msdanalyzer object
            temp_simSPT_track_msd = msdanalyzer(2,'um','seconds');
                    
            for iter = 1:length(now_track_cell)
                now_track = now_track_cell{iter};                        
                temp_simSPT_track_msd = temp_simSPT_track_msd.addAll({[now_track.t,now_track.x,now_track.y]});
            end
            
            % calculate the MSD of all imported tracks
            temp_simSPT_track_msd = temp_simSPT_track_msd.computeMSD();
            
            % Append mean MSD into the table mean_T
            for iter = 1:length(now_track_cell)
                % Write area into a table
                now_track = now_track_cell{iter}; 
                new_line = table(pix_filename(iter),now_track.trajecotry_newIdx(1),target_exposure(idT),simSPT_D_num(iSPT),...
                    temp_simSPT_track_msd.msd{iter}(2,2),... % mean MSD at 1*dT
                    mask_merge_z_area(iter),mask_pix_merge_z_area(iter),...
                    ...fitresult(iter,8),fitresult(iter,9),...
                    JumpStartToEnd(iter));
                writetable(new_line,fullfile(SummaryTable_savePath,sprintf('NoNoise_SmoothArea_PixelArea.csv')),'WriteMode','append',...
                            'WriteVariableNames',false,'WriteRowNames',true); 
            end
        else
            % Append mean MSD into the table mean_T
            for iter = 1:length(now_track_cell)
                % Write area into a table
                now_track = now_track_cell{iter}; 
                new_line = table(pix_filename(iter),now_track.trajecotry_newIdx(1),target_exposure(idT),simSPT_D_num(iSPT),...
                    ...temp_simSPT_track_msd.msd{iter}(2,2),... % mean MSD at 1*dT
                    mask_merge_z_area(iter),mask_pix_merge_z_area(iter),...
                    ...fitresult(iter,8),fitresult(iter,9),...
                    JumpStartToEnd(iter));
                writetable(new_line,fullfile(SummaryTable_savePath,sprintf('NoNoise_SmoothArea_PixelArea.csv')),'WriteMode','append',...
                            'WriteVariableNames',false,'WriteRowNames',true); 
            end
        end

        % save log 
        fprintf(fileID,'%3s\t%3d\t%12d\n',simSPT_D(iSPT),target_exposure(idT),length(now_track_cell));
    
    end   
end

fclose(fileID);
delete(gcp('nocreate')); % delete the current parrallel pool
fprintf('============================ All completed! ============================\n');

sourceFile = mfilename('fullpath'); % Define the source and target file paths
copyfile([sourceFile '.m'], img_pair_saveDir);

% %% AUXILIARY
% function PSF_BoundxyI = make_BoundxyI(Xgrid,Ygrid,merge_z)
%     tempx = Xgrid([1,end],2:31); % CHANGE HERE if pixelized image dimension changed 
%     tempy = Ygrid([1,end],2:31); % CHANGE HERE if pixelized image dimension changed 
%     tempz = merge_z([1,end],2:31); % CHANGE HERE if pixelized image dimension changed 
%     PSF_BoundxyI = [Xgrid(:,1), Ygrid(:,1), merge_z(:,1); tempx(:),tempy(:),tempz(:);Xgrid(:,end), Ygrid(:,end), merge_z(:,end)];
% end
