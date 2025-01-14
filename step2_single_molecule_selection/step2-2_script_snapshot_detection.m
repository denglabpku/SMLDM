% This is a script to do batch-mode PSF extraction
% Input ND2 requires MPALM-BULK-MPALM... image sequence, and should be properly arranged
% see 2024/20240125_U2OS_HaloRPB1_HMSiR_PA646/HMSiR10nM_TMR/Clust02 for example
% currently do not support multiple cells in one FOV

%% DESCRIPTION
% Use bfmatlab version >6.6.0 to directly load nd2 file without error


%% Zuhui Wang
%% 2023/02/20

%% WARNING NOTE
% Warning: Duplicate data points have been detected and
% removed. can be safely ignored.
% NOTE: FINAL SNRI IS NOT UPDATED AFTER ADDMTT

%% INTENTIONALLY KEEP BLANK
%% INTENTIONALLY KEEP BLANK
%% INTENTIONALLY KEEP BLANK
clc; clear; close all;
addpath(genpath('/home/wdeng_pkuhpc/lustre3/wangzuhui/apps/matlab/psfanalysis')) % add bio-formats MATLAB Toolbox to the search path
%% DEFINE PARAMETERS

% >>>>>>>>>>>>>>>>>>>> NUCLEUS SELECTION (LARGE FOV ONLY) >>>>>>>>>>>>>>> %
has_bulk_imgs = false;
skip_drawROI = false; % set true if already done
draw_fullROI = true;
MPALM_CHANNEL = 1; % channel of mPALM. 1 for left cam, 2 for back cam.
WideField_subfolder = 'BF';
% <<<<<<<<<<<<<<<<<<<< NUCLEUS SELECTION (LARGE FOV ONLY) <<<<<<<<<<<<<<< %

% >>>>>>>>>>>>>>>>>>>> MOTION BLUR DETECTION PARAMETERS >>>>>>>>>>>>>>>>>>>>> %
impars.PixelSize = 110; % nm per pixel
impars.EmissionWavelength = 664; % wavelength in nm; consider emission max and filter cutoff; 580 for mEos3.2 & mEosEM; 582 for Blinking TMR;664 for PA-JF646; 594 for BD566blink; 668 for HaloX-SiR
impars.PhotonConversionRate = 0.9398923112; % (0.5e/ADU, EMGain300, QE at 660nm:0.9398923112) (0.5e/ADU, EMGain300, QE at 582nm:0.9696909185, QE at 594nm:0.96593569835; QE at 668nm: 0.934925852)
impars.ExposureTime = 30.5; % in milliseconds

% General parameters
filter.MaxTotalPhoton = Inf; % Total photon number exceed this value will be discarded, 2000 is a good value to remove out-of-focus blurry signal; Use ThunderSTORM to keep only one loc MB, so not necessary.
filter.UNet_model_tif = 'UNet_mask_MBX_20240620_2035_epoch20_Ch1.tif';
filter.ThunderSTORM_saveSubfolder = 'ThunderSTORM';
filter.useLocFilter = true; % Use ThunderSTORM localization to filter molecule signals
% <<<<<<<<<<<<<<<<<<<< MOTION BLUR DETECTION PARAMETERS <<<<<<<<<<<<<<<<<<<<< %

% DEFINE INPUT AND OUTPUT
% where you save ND2 file
input_path = '/home/wdeng_pkuhpc/lustre3/wangzuhui/raw_data/YinChao/20241127_PA646_TMR_CJ9_HomoKI_Halo-RPB1_MPALM/';

% the parent folder of motion blur detection and analysis
output_path = '/home/wdeng_pkuhpc/lustre3/wangzuhui/processed_data/Yinchao/20241127_PA646_TMR_CJ9_HomoKI_Halo-RPB1_MPALM/';

% ND file name
input_rawND2_prefix = {...
    '20241127_Cell05_CJ9_Halo_RPB1',
    '20241127_Cell06_notreat_CJ9_Halo_RPB1',
    };

% folder path of where to save main results
ROI_infoPath = cellfun(@(x) fullfile(output_path,'roi_files',[x '_roi_metadata.mat']),input_rawND2_prefix,'UniformOutput',false);
mkdir(fullfile(output_path,'roi_files'));

%% Nuclear ROI selection
for iSamp = 1:length(input_rawND2_prefix)

    if has_bulk_imgs
        WideField = dir(fullfile(input_path,WideField_subfolder,sprintf('Clust%02d',input_ND2_clustIdx(iSamp)),'*.nd2'));
    end
        
    %% SEGMENT NUCLEUS MANNUALLY AND EXPORT THE NUCLEUS ROI INFO
    % [Bug] if unexpected exit before finish, nucleus_good_roi.csv can not be
    % written.
    if ~skip_drawROI
        if draw_fullROI
            disp('... No need for ROI drawing, generate a fake ROI cover whole FOV ...'); 
            
            input_rawND2_list = dir(fullfile(input_path,'*.nd2'));
            Filenames = {input_rawND2_list.name};
            valid_idx = contains(Filenames,input_rawND2_prefix(iSamp)); 
            input_rawND2_list = input_rawND2_list(valid_idx);            
            [currentImage,~] = MemoryEfficientND2reader_oneFrame(fullfile(input_rawND2_list(1).folder, input_rawND2_list(1).name),1);
            [ImHeight,ImWidth] = size(currentImage); clear currentImage;
            % img_stack_cell_array = bfopen(fullfile(WideField(1).folder,WideField(1).name));   
            % currentImage_Widefield = img_stack_cell_array{1}{1,1}; % 3d image matrix
            % [ImHeight,ImWidth] = size(currentImage_Widefield); clear currentImage_Widefield;
            roi_info_nuc = {};
            % for iBulk = 1:length(WideField)
            for iBulk = 1:length(input_rawND2_list)            
                roi_info_nuc(iBulk) = {[0,0; ImWidth,0; ImWidth, ImHeight; 0,ImHeight]};
            end
            save(ROI_infoPath{iSamp},'roi_info_nuc');

        else

            previousROI = []; % Initializing the previous ROI
            roi_info_nuc = {};

            for iBulk = 1:length(WideField)
                close all
                img_stack_cell_array = bfopen(fullfile(WideField(iBulk).folder,WideField(iBulk).name));   
                currentImage_Widefield = img_stack_cell_array{1}{1,1}; % 3d image matrix

                % Display the image
                fig_ROI = figure('Position',[680 318 556 505],'Name',sprintf('Bulk image: #%02d', iBulk));
                imshow(currentImage_Widefield,[],'Border','tight','InitialMagnification',1000);
                
                ax_ROI = fig_ROI.CurrentAxes;

                if iBulk == 1
                    hI = imcontrast;
                    % Let the user draw the ROI
                    hROI = drawpolygon(ax_ROI,'LabelVisible','on',...
                                    'FaceAlpha', 0,'LineWidth', 2,'Color','green');
                    wait(hROI);
                    position = hROI.Position; % Each row represents the [x y] coordinates of a ROI vertex
                    fig_CLim = [ax_ROI.CLim(1) ax_ROI.CLim(2)];
                else
                    % Reuse last ROI
                    ax_ROI.CLim = fig_CLim;
                    % hI = imcontrast;
                    previousROI = roi_info_nuc{iBulk-1};
                    hROI = drawpolygon(ax_ROI,'Position',previousROI,'LabelVisible','on',...
                                    'FaceAlpha', 0,'LineWidth', 2,'Color','green',...
                                    'InteractionsAllowed','translate');
                    f = msgbox("Finish ROI adjustment?"); f.Position = [921.7500 406.5000 150 51.7500];
                    uiwait(f);
                    position = hROI.Position;
                end
                
                roi_info_nuc(iBulk) = {position};
            end

            save(ROI_infoPath{iSamp},'roi_info_nuc');
        end
        
        fprintf('All ROI selection finished!\n');
        % fprintf('Nucleus ROI with good evenness saved to : %s\n',fullfile(output_path,'roi_files','nucleus_good_roi.csv'));
        close all;
        
    end % IF skip_drawROI
end

%% INTENTIONALLY KEEP BLANK
%% INTENTIONALLY KEEP BLANK

%% MOTION BLUR DETECTION

disp('============================================')
disp('            _   _ __  __ ______  __')
disp('           | | | |  \/  | __ ) \/ /')
disp('           | | | | |\/| |  _ \\  /')
disp('           | |_| | |  | | |_) /  \')
disp('            \___/|_|  |_|____/_/\_\')
disp('============================================')
disp('Motion Blur Extraction from Single Molecules');
disp(' Designed by Zuhui Wang, Peking University');
disp('     Current version: V2p2 MPALM-BULK sequence');
disp('============================================')
disp('============================================')

NumWorkers = 28;
parpool('local', NumWorkers)

for iSamp = 1:length(input_rawND2_prefix)
    tic;

    fprintf('Start processing the current image sequence: %s\n',input_rawND2_prefix{iSamp});    
    clearvars roi_info_nuc imgs_2d_matrix ImHeight ImWidth tot_img_num mask_img loc_data tot_cell_num

    input_rawND2_list = dir(fullfile(input_path,'*.nd2'));
    Filenames = {input_rawND2_list.name};    
    valid_idx = contains(Filenames,input_rawND2_prefix(iSamp)); 
    input_rawND2_list = input_rawND2_list(valid_idx);

    UNet_model_folder = cellfun(@(x) fullfile(output_path,x(1:end-4)),{input_rawND2_list.name},'UniformOutput',false);
    tab_psf_fitresult = table();    

    % create subfolder name for each cell
    clearvars cell_save_folder
    cell_save_folder = fullfile(output_path,sprintf('%s_Cell%02d',input_rawND2_prefix{iSamp},1));
    if ~exist(cell_save_folder,'dir')
        mkdir(cell_save_folder)
    end

    %============================ Read ROI =================================%
    load(ROI_infoPath{iSamp},'roi_info_nuc')
    tot_slice_num = length(input_rawND2_list);
    impars.tot_img_num = 0;

    for iSlice = 1:tot_slice_num
        fprintf('Processing the #%d/#%d slice in this sequence...\n',iSlice,tot_slice_num);

        cell_PSF.xyI = [];
        cell_PSF.BoundxyI = [];
        cell_PSF.label = []; % can use to draw masks
        cell_PSF.boundLabel = []; % can use to draw mask outline
        
        %============================ Read ND2 Image =============================%
        fprintf('Reading ND2 image ...\n');    
        temp_image_file = fullfile(input_rawND2_list(iSlice).folder, input_rawND2_list(iSlice).name);
        [imgs_2d_1stframe,tot_img_num] = MemoryEfficientND2reader_oneFrame( temp_image_file,  1);
        impars.tot_img_num = impars.tot_img_num + tot_img_num;        
        [ImHeight,ImWidth] = size(imgs_2d_1stframe);
        impars.ImHeight = ImHeight;
        impars.ImWidth = ImWidth;        

        %====================== READING LOCALIZATION FROM THUNDERSTORM ======================%
        disp('... Reading ThunderStorm localization');
                
        loc_data_table = readtable(fullfile(output_path,filter.ThunderSTORM_saveSubfolder,sprintf('%s_C%d.csv',input_rawND2_list(iSlice).name(1:end-4),MPALM_CHANNEL)),'VariableNamingRule','preserve');
        loc_data_table.xPix = loc_data_table.("x [nm]")/impars.PixelSize+1;
        loc_data_table.yPix = loc_data_table.("y [nm]")/impars.PixelSize+1;
            
        %==================== Keep MTT locs inside ROI =======================%
        fprintf('... Apply nucleus mask\n');
        nuclear_mask = logical(poly2mask(roi_info_nuc{iSlice}(:,1),roi_info_nuc{iSlice}(:,2),ImHeight,ImWidth));
        
        yPix_rounded = round(loc_data_table.yPix);
        xPix_rounded = round(loc_data_table.xPix);
        linear_indices = sub2ind(size(nuclear_mask), yPix_rounded, xPix_rounded);
        in_roi_array = nuclear_mask(linear_indices);
        loc_data_nuclear_table = loc_data_table(logical(in_roi_array),:);
        par_locTable = cell(tot_img_num,1);
        for i_plane = 1: tot_img_num
            par_locTable{i_plane} = loc_data_nuclear_table(loc_data_nuclear_table.frame == i_plane,:);
        end


        % Prepare parfor temporary variable in each iteration
        cell_PSF_xyI = cell(tot_img_num,1); % pixel coordinates and the corresponding intensity of PSF masks
        cell_PSF_BoundxyI = cell(tot_img_num,1); % pixel coordinates and the corresponding intensity of PSF mask boundaries.
        cell_PSF_SNRI = cell(tot_img_num,1); 
        cell_PSF_label = cell(tot_img_num,1);
        cell_PSF_boundLabel = cell(tot_img_num,1);
        cell_psftable = cell(tot_img_num, 1);

        temp_temp_UNet_model_folder = UNet_model_folder{iSlice};        
        temp_UNet_model_tif = filter.UNet_model_tif;

        parfor i_plane = 1:tot_img_num
            %============================ Read ND2 Image =============================%
            [imgs_2d_matrix,~] = MemoryEfficientND2reader_oneFrame(temp_image_file, i_plane);
            imgs_2d_matrix = double(imgs_2d_matrix);
        
            %======================= Read UNet generated masks =======================%
            % fprintf('Reading UNet generated masks ...\n');
            % Initialize an empty 3D array to hold the images
            mask_img = imread(fullfile(temp_temp_UNet_model_folder,temp_UNet_model_tif),i_plane);
                            
            [temp_cell_PSF_xyI,temp_cell_PSF_BoundxyI,temp_cell_PSF_SNRI,PSF_label,PSF_boundLabel] = MB_Extraction_inROI_ThunderStorm(...
            imgs_2d_matrix,mask_img,nuclear_mask,filter,par_locTable{i_plane},impars);

            % par_cell_ToCloseLoc{i_plane} = temp_par_ToCloseLoc;
            cell_PSF_xyI{i_plane} = temp_cell_PSF_xyI;
            cell_PSF_SNRI{i_plane} = temp_cell_PSF_SNRI;
            cell_PSF_BoundxyI{i_plane} = temp_cell_PSF_BoundxyI;
            cell_PSF_label{i_plane} = PSF_label;
            cell_PSF_boundLabel{i_plane} = PSF_boundLabel;
            
            %=============================== Ellipticity fitting ===============================%   
            [fitresult, ~, ~]= PARA_GAUSS2DFIT(vertcat(temp_cell_PSF_xyI),vertcat(temp_cell_PSF_BoundxyI));
            fitresult(:,8) = ones(height(fitresult),1)*i_plane; % indicate frame number
            fitresult(:,9) = (1:height(fitresult))'; % indicate PSF number in this frame

            % append ellipticity
            fitresult(:,10) = max(fitresult(:,3:4),[],2)./min(fitresult(:,3:4),[],2);
            %fitresult = [amp, ang, sx, sy, xo, yo, zo, frame#, PSFLabel#, ellipticity]
            
            temp_SNR = vertcat(temp_cell_PSF_SNRI);
            temp_COV = cellfun(@(x) std(x(:,3))/mean(x(:,3)),vertcat(temp_cell_PSF_xyI));
            temp_area = cellfun(@height,vertcat(temp_cell_PSF_xyI));
            temp_totalphoton = cellfun(@(x)sum((x(:,3)-500)/(0.5*300*impars.PhotonConversionRate)),vertcat(temp_cell_PSF_xyI));
            temp_angle = fitresult(:,2);
            for i = 1:height(fitresult) % Correct the fitted angle to make the angle conter-clockwise starting from x axis.
                if fitresult(i,3) <= fitresult(i,4) % check which is the primary axis of the fitting angle
                    if 90-fitresult(i,2) <0 
                        temp_angle(i) = 90-fitresult(i,2)+180;
                    else
                        temp_angle(i) = 90-fitresult(i,2);
                    end
                else
                    temp_angle(i) = 180-fitresult(i,2);
                end
            end
            
            %============================ PREPARE EXPORT FINAL RESULT ============================== %
            temp_tab_psf_fitresult = table(fitresult(:,8)+(tot_img_num*(iSlice-1)),fitresult(:,9),fitresult(:,5),fitresult(:,6),temp_totalphoton,fitresult(:,1),fitresult(:,7),100*log10(fitresult(:,10)),temp_SNR,temp_COV,temp_area,temp_angle,...
                'VariableNames',{'Frame','PSF_idx','Xpos','Ypos','TotalPhoton','Intensity','Background','EllipticityIdx','SNR','COV','Area','Angle'});                  
            
            cell_psftable{i_plane} = temp_tab_psf_fitresult;            
        end

        cell_PSF.xyI = cell_PSF_xyI;
        cell_PSF.BoundxyI = cell_PSF_BoundxyI;
        cell_PSF.label = cell_PSF_label; % can use to draw masks
        cell_PSF.boundLabel = cell_PSF_boundLabel; % can use to draw mask outline.

        % Concatenate results after the parfor loop
        temp_tab_psf_fitresult = vertcat(cell_psftable{:});
        temp_tab_psf_fitresult = sortrows(temp_tab_psf_fitresult, {'Frame','PSF_idx'}); % Sort rows for PSF_idx first, then for Frame
        tab_psf_fitresult = [tab_psf_fitresult;temp_tab_psf_fitresult];  

        % Save workspace if everything looks good
        save(fullfile(cell_save_folder,sprintf('Blurdata_%s_Slice%02d.mat',filter.UNet_model_tif(1:end-4),iSlice)),...
            'input_path', 'output_path','input_rawND2_prefix','UNet_model_folder','filter','iSlice',...
            'impars',...
            'cell_PSF',...
            '-v7.3');   
        
        clearvars cell_PSF cell_PSF_xyI cell_PSF_BoundxyI cell_PSF_label cell_PSF_boundLabel cell_PSF_SNRI
    end

    writetable(tab_psf_fitresult,fullfile(cell_save_folder,sprintf('Fitresult_%s.csv',filter.UNet_model_tif(1:end-4))))
    
    clearvars tab_psf_fitresult
    toc;

end
delete(gcp('nocreate')); % delete the current parrallel pool
fprintf('============================ All completed! ============================\n');

