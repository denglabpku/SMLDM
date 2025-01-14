function [temp_cell_PSF_xyI,temp_cell_PSF_BoundxyI,temp_cell_PSF_SNRI,PSF_label,PSF_boundLabel] = MB_Extraction_inROI_ThunderStorm(...
    I,raw_mask,nuclear_mask,filter,par_locTable,impars)
%MB_Extraction Extract motion-blur(called PSF here) from UNet mask
%   Detailed explanation goes here    
%% Zuhui Wang
%% 2023/06/11

%Initiate some variables
temp_cell_PSF_xyI = cell(0,1);
temp_cell_PSF_BoundxyI = cell(0,1);
temp_cell_PSF_SNRI = [];
PSF_boundLabel = zeros(size(I));

%Fill holes
BWdfill = imfill(raw_mask, 'holes');

%Clear Border Objects
BWnoboard = imclearborder(BWdfill,4);

% Below option are better not used.U-Net is good enough 
% %Option1 Unbridge the Objects if only one pixel connected, change bwlabel
% %connectivity to 4
% BWnoboard = ~bwmorph(~BWnoboard,'bridge');
% %Option2 Watershed
% L = watershed(-bwdist(~BWnoboard));
% BWnoboard = BWnoboard&L~=0;
% %Option3
% Post-processing with Morphological Operations (better than using disk structure element)
% se90 = strel('line',2,90); %Define vertical structuring elements
% se0 = strel('line',2,0); %Define horizontal structuring elements
% J1 = imopen(BWnoboard,se90);
% BWnoboard = imopen(J1,se0);

%Label connected components in 2-D binary image; use 8 for watershed otherwise use 4
PSF_label = bwlabel(BWnoboard, 4); 

%% PSF QUALITY CHECK
%============ Size and ROI =======================
% label PSF as background if its size is not corrrect
% Set the threshold for object sizes
% PSF size filter, Remove PSFs that size are abnormally large
% (usually due to close adjacent signals) or small (noise)
% This is the most efficient way to remove the rest abnormally
% large PSFs, hard to eliminate thru other filters or methods.

% Get the region properties of the labeled regions
stats = regionprops(PSF_label, ["Area","Centroid"]);
for i = 1:numel(stats)
    % keep label only inside ROI
    round_xy = round(stats(i).Centroid);
    if nuclear_mask(round_xy(2),round_xy(1)) == 0
        PSF_label(PSF_label == i) = 0;
    end
%     % Loop through each region and fill too small or too big ones with zeros
%     if stats(i).Area < filter.PSFSz_bottom || stats(i).Area > filter.PSFSz_top
%         PSF_label(PSF_label == i) = 0;
%     end
    % % Loop through each region and fill too small ones with zeros
    % if stats(i).Area < filter.PSFSz_bottom 
    %     PSF_label(PSF_label == i) = 0;
    % end
    
    % discard PSF with too high total photon number
    % EXTRACT PSF REGION
    [PSF_rIdx, PSF_cIdx] = find(PSF_label == i);
    PSF_lIdx = sub2ind(size(I),PSF_rIdx,PSF_cIdx); % Convert subscripts to linear indices       
    PSF_Intensity = double(I(PSF_lIdx));
    photon = (PSF_Intensity-500)/(0.5*300*impars.PhotonConversionRate);
    if sum(photon) > filter.MaxTotalPhoton
        PSF_label(PSF_label == i) = 0;
    end
    clearvars PSF_rIdx PSF_cIdx PSF_lIdx;
end
% Relabel the image starting from 1
PSF_label = bwlabel(PSF_label,4);


% %============ ThunderStorm locs number in one PSF =======================

PSF_num = max(PSF_label(:)); 
if PSF_num ~= 0
    % Loop over all points and label them based on their location in the mask image
    points = [par_locTable.xPix,par_locTable.yPix]; % [x y] starting from 1, in matlab image coordinate
    ThunderStormLoc_inPSF_labels = zeros(size(points, 1), 1);
    for i = 1:height(points)
        x = round(points(i, 1)); % assume points are stored as (x,y) pairs
        y = round(points(i, 2));
        ThunderStormLoc_inPSF_labels(i) = PSF_label(y, x);
        % if mask_binary(y, x) == 1 % point belongs to a mask region
        %     ThunderStormLoc_inPSF_labels(i) = PSF_label(y, x);
        % else % point does not belong to any mask region
        %     ThunderStormLoc_inPSF_labels(i) = 0;
        % end
    end
    
    for PSF_idx = 1:PSF_num
        % Current PSF contains MTT locs number
        ThunderStormLoc_inOnePSF = sum(ThunderStormLoc_inPSF_labels == PSF_idx);
        % Remove PSF that contains 2 MTT locs
        % if ThunderStormLoc_inOnePSF > 1
        if ThunderStormLoc_inOnePSF ~= 1
            PSF_label(PSF_label == PSF_idx) = 0;

            % % label MTT locs in same PSF
            % in_same_PSF = ismember(PSF_idx, ThunderStormLoc_inPSF_labels);

            % temp_par_ToCloseLoc(in_same_PSF) = true; % A logic to tell if both MTT are too close, used in PSF_ADDMTT
            % temp_cell_PSF_MTTxy(PSF_idx,:) = [0,0]; % No MTT locs inside the detected PSF

        % elseif ThunderStormLoc_inOnePSF == 1
        %     temp_cell_PSF_MTTxy(PSF_idx,:) = [par_ctrsX(loc_idx_range(ThunderStormLoc_inOnePSF))+1, par_ctrsY(loc_idx_range(ThunderStormLoc_inOnePSF))+1];
        % elseif ThunderStormLoc_inOnePSF == 0
        %     temp_cell_PSF_MTTxy(PSF_idx,:) = [0,0]; % No MTT locs inside the detected PSF
        end
    end
end

% Relabel the image starting from 1
PSF_label = bwlabel(PSF_label,4);

% %============ (To-do) Keep PSF inside ROI =======================


%% SAVE PSF
%================== Begin saving every good PSF ==========================

PSF_num = max(PSF_label(:)); 
if PSF_num ~= 0
    % Save each PSF into cell array
    PSF_rcIMat = cell(PSF_num,1);%To save PSF region pixels
    PSF_rcIBMat = cell(PSF_num,1); %To save PSF boundary region pixels
    
    % Get boundary of each PSF
    opt_boundaryMask = 2; %1 to use boundarymask(); 2 to use bwperim(). Don't change.
    if opt_boundaryMask == 1
        %Option1: Use boundarymask(), wider boundary width, but when two PSF
        %are too close, boundary may merge
        PSF_maskBound = boundarymask(BWnoboard,4);%Find boundary pixels of the masks
        PSF_labelTemp = PSF_label;
        PSF_labelTemp(PSF_maskBound~=1) = 0;
        PSF_boundLabel = PSF_labelTemp;
    elseif opt_boundaryMask == 2
        %Option2: Use bwperim(), 1 pixel boundary width, but can separate close
        %PSF well 
        PSF_maskBound = bwperim(PSF_label,4); % use 4 for later visualization purpose only. find the perimeter pixels of objects in the input image BW, so still part of BW
        PSF_labelTemp = PSF_label;
        PSF_labelTemp(PSF_maskBound~=1) = 0;
        PSF_boundLabel = PSF_labelTemp;   
    end
    
    for PSF_idx = 1:PSF_num

        % EXTRACT PSF REGION
        [PSF_rIdx, PSF_cIdx] = find(PSF_label == PSF_idx);
        PSF_lIdx = sub2ind(size(I),PSF_rIdx,PSF_cIdx); % Convert subscripts to linear indices       

        %Enumerate the pixel coordinate and its pixel intensity in
        %left-to-rght and column/x-coord first way (see bwlabel for detail)
        %Store each PSF coordinate and intensity in cell array
        PSF_rcIMat{PSF_idx} = [double(PSF_cIdx) double(PSF_rIdx) double(I(PSF_lIdx))];

        clearvars PSF_rIdx PSF_cIdx PSF_lIdx; %Recover valuable variable

        % EXTRACT EXACT GRADIENT MASK PSF REGION BOUNDARY
        [PSF_rIdx, PSF_cIdx] = find(PSF_boundLabel == PSF_idx);
        PSF_lIdx = sub2ind(size(I),PSF_rIdx,PSF_cIdx);
        PSF_rcIBMat{PSF_idx} = [double(PSF_cIdx) double(PSF_rIdx) double(I(PSF_lIdx))];
        
        clearvars PSF_rIdx PSF_cIdx PSF_lIdx;

    end
    
    %%% Option3: Calculate SNR of each PSF using noise from whole image
    % SNR definition from Serge et al, Nat Method, 2008.    
    PSF_maxI = cellfun(@(x)max(x(:,3)),PSF_rcIMat); %Max of x, y coord and intensity of PSF
%     PSF_maxI = cell2mat(PSF_maxI);
    % variance should be the noise variance, do not include the
    % signal!!!!!! Pixel Intensity = signal + noise + offset (Serge, 2008 MTT)
    % PSF_SNRI = 10.*log10((PSF_maxI(:,3).^2)./PSF_varianceBoundI(:,3)); %SNR in decibels
    % PSF_SNRI = 10.*log10(((PSF_maxI-median(I,'all')).^2)./filter.noise_std^2); %SNR in decibels
    
    % Use mean-threshold to dynamically estimate background
    background_pix = I(I <= mean(I,'all'));
    background = mean(background_pix,"all");
    noise = std(background_pix,0,'all');
    PSF_SNRI = 10.*log10((PSF_maxI-background).^2./noise^2); %SNR in decibels
    

    temp_cell_PSF_xyI = PSF_rcIMat;
    temp_cell_PSF_BoundxyI = PSF_rcIBMat;
    temp_cell_PSF_SNRI = PSF_SNRI;
end

end