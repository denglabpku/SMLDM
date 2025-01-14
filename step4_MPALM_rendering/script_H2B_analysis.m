% This script is to analyze H2B MPALM

clear;close all;clc;

% Load dataset
SRFile.folder = 'F:\PROCESS-SPT\2023\20231013_U2OS_H2BOnly_HP1-H2B_mPALM\EF1a-H2B\20231013_Clust02_EF1a-H2B_642_25uMPA646_30p5ms_20kframe_02_Cell01';
SRFile.name = 'drift_corrected_Ch1_SR_pred_v3.csv';

% SRFile.folder = 'F:\PROCESS-SPT\2023\20231013_U2OS_H2BOnly_HP1-H2B_mPALM\EF1a-H2B\20231013_Clust01_EF1a-H2B_642_25uMPA646_30p5ms_20kframe_01_Cell01';
% SRFile.name = 'drift_corrected_SR_pred_v3.csv';

tab_SR_pred = readtable(fullfile(SRFile.folder,SRFile.name)); % Xpos, Ypos unit already nm if use drift-corrected table

Blurdata_mat = load(fullfile(SRFile.folder, 'Blurdata_UNet_mask_MBX_20240620_2035_epoch20_Ch1.mat'), 'impars');
params.ImWidth = Blurdata_mat.impars.ImWidth;
params.ImHeight = Blurdata_mat.impars.ImHeight;

params.dxy = 30; % rendered pixel size
params.sigma_render = 50;
params.camPixelSize_Ch1 = 110; % camera pixel size, unit nm.
params.Ch1_loc_XYLim = [0,params.camPixelSize_Ch1*params.ImWidth,0,params.camPixelSize_Ch1*params.ImHeight];


%% Modified HMRFsegmentation function to segment chromatin regions
dxy = params.dxy; % rendered pixel size
sigma_render = params.sigma_render;

Xmin = params.Ch1_loc_XYLim(1);Xmax=params.Ch1_loc_XYLim(2);Ymin = params.Ch1_loc_XYLim(3);Ymax=params.Ch1_loc_XYLim(4);
Edges{1}=Xmin:dxy:Xmax+dxy; % make edge cover all data points
Edges{2}=Ymin:dxy:Ymax+dxy;


% Xmin = min(tab_SR_pred.Xpos);Xmax=max(tab_SR_pred.Xpos);Ymin = min(tab_SR_pred.Ypos);Ymax=max(tab_SR_pred.Ypos);
% Edges{1}=Xmin:dxy:Xmax+dxy; % make edge cover all data points
% Edges{2}=Ymin:dxy:Ymax+dxy;

% Create mask cover the nucleus
% bin every localizations into grid, similar to Im = hist3([Xpos',Ypos'],'Edges',Edges);
[Im,Xedges,Yedges,binX,binY] = histcounts2(tab_SR_pred.Xpos,tab_SR_pred.Ypos,Edges{1},Edges{2}); 
TempX=-round(3*sigma_render/dxy)*dxy:dxy:round(3*sigma_render/dxy)*dxy;
TempY=-round(3*sigma_render/dxy)*dxy:dxy:round(3*sigma_render/dxy)*dxy;

ConVecX = exp(-0.5*(TempX/sigma_render).^2);
ConVecX=ConVecX/sum(ConVecX);
ConVecY = exp(-0.5*(TempY/sigma_render).^2);
ConVecY=ConVecY/sum(ConVecY);

Im2 = conv2(ConVecX,ConVecY,Im,'same');
ImPALM = transpose(Im2/dxy/dxy*10^6);

figure;
imshow(ImPALM,[],color=hot);hIc = imcontrast(); uiwait(hIc);
% % Draw a polygon ROI
% hROI = drawpolygon('LineWidth', 2, 'Color', 'r');
% % Get the X and Y coordinates of the vertices of the polygon
% roiVertices = hROI.Position;
% 
% figure;
% imshow(transpose(Im),[],color=hot);hIc = imcontrast(); uiwait(hIc);
% % Draw a polygon ROI
% hROI = drawpolygon('Position',roiVertices,'LineWidth', 2, 'Color', 'r');


% HMRF segmentation (currently for H2B only) using immobile molecule
immobile_table = tab_SR_pred(tab_SR_pred.D <= 0.05,:);
% bin every localizations into grid, similar to Im = hist3([Xpos',Ypos'],'Edges',Edges);
[immobileIm,~,~,~,~] = histcounts2(immobile_table.Xpos,immobile_table.Ypos,Edges{1},Edges{2}); 
TempX=-round(3*sigma_render/dxy)*dxy:dxy:round(3*sigma_render/dxy)*dxy;
TempY=-round(3*sigma_render/dxy)*dxy:dxy:round(3*sigma_render/dxy)*dxy;

ConVecX = exp(-0.5*(TempX/sigma_render).^2);
ConVecX=ConVecX/sum(ConVecX);
ConVecY = exp(-0.5*(TempY/sigma_render).^2);
ConVecY=ConVecY/sum(ConVecY);

immobileIm2 = conv2(ConVecX,ConVecY,immobileIm,'same');

immobileImPALM = transpose(immobileIm2/dxy/dxy*10^6);

nclust = 7; % number of total chromatin segmentation, default 7.
xylim = [Xmin, Xmax, Ymin, Ymax];

[HMRFseg, setting] = HMRFsegmentation(immobile_table,xylim,[],nclust,dxy,sigma_render,0.1,immobileImPALM);

% merge 7 class into 2 class (non-CD and CD)
HMRFseg.simple_img_class = HMRFseg.img_class;

HMRFseg.simple_img_class(ismember(HMRFseg.img_class,[1:3])) = 0;
HMRFseg.simple_img_class(ismember(HMRFseg.img_class,[4:7])) = 1;

full_partition = reshape(HMRFseg.img_class, size(HMRFseg.nucleus_mask)); full_nucleus_partition = transpose(full_partition);

simple_nucleus_partition = reshape(HMRFseg.simple_img_class, size(HMRFseg.nucleus_mask)); simple_nucleus_partition = simple_nucleus_partition';

cmap_all = uint8([4, 0, 0; 56, 46, 142; 137, 48, 141; 215, 31, 40; 239, 127, 25; 244, 191, 27; 244, 237, 70; 255, 255, 255]);
cmap_CD = uint8([4, 0, 0; 56, 46, 142; 215, 31, 40;255, 255, 255]);

% make this as a montage
figure('Units','normalized','Position',[0.1 0.1 0.5 0.5]);
tiledlayout(1,2,'TileSpacing','none','Padding','tight');

ax1 = nexttile; hold on;
imagesc(ax1,full_nucleus_partition)
% hROI = drawpolygon('Position',roiVertices,'LineWidth', 2, 'Color', 'r');
colormap(ax1,cmap_all)
axis(ax1,'off')
axis image ij
daspect(ax1,[1, 1, 1]);
title(ax1,'HMRF chromatin segmentation')

ax2 = nexttile;
imagesc(ax2,simple_nucleus_partition)
colormap(ax2,cmap_CD)
axis(ax2,'off')
daspect(ax2,[1, 1, 1]);
title(ax2,'Binary CD and non-CD')

%% Plot

% Average D at each chromatin class
% Initialize accumulators for sum of D values and counts per label
sumD = zeros(nclust, 1);
countD = zeros(nclust, 1);
sumSquaredD = zeros(nclust, 1);

sumD_slow = zeros(nclust, 1);
countD_slow = zeros(nclust, 1);

D_values = cell(nclust, 1); % Ensure this matches the number of clusters
% Iterate over each data point
for i = 1:numel(binX)
    if binX(i) > 0 && binY(i) > 0 && full_partition(binX(i), binY(i)) > 0
        % Get the label from full_partition
        label = full_partition(binX(i), binY(i));
               
        % Accumulate D and count for the current label
        dValue = tab_SR_pred.D(i);
        sumD(label) = sumD(label) + dValue;
        sumSquaredD(label) = sumSquaredD(label) + dValue^2;
        countD(label) = countD(label) + 1;
        
        % Store D value in a cell array
        D_values{label} = [dValue, D_values{label}];
        

        if tab_SR_pred.D(i) < 0.05
            sumD_slow(label) = sumD_slow(label) + dValue;
            countD_slow(label) = countD_slow(label) + 1;
        end
    end
end


% Calculate average D for each label
allLoc_averageD = sumD ./ countD;
slowLoc_averageD = sumD_slow./countD_slow;

% Detection fraction
allLoc_Fcount = countD/sum(countD);

% Calculate standard deviation for each label
stdDevD = sqrt((sumSquaredD ./ countD) - (averageD.^2));
figure;
errorbar(allLoc_averageD,stdDevD)
axis padded

J1 = customcolormap_preset('pasteljet');
J2 = customcolormap_preset('red-white-blue');

figure; tiledlayout(1,2);
nexttile;heatmap(allLoc_averageD,'Colormap',J1);title('Global average D');
nexttile;heatmap(slowLoc_averageD,'Colormap',J1);title('Immobile average D');

figure;
% histogram(ImPALM(:),'Normalization','percentage')
heatmap(allLoc_Fcount,'Colormap',J2);title('Global count fraction')


% % --------Intensity distribution and Gaussian mixture model--------
% fig1 = figure;
% hold on
% h = histogram(HMRFseg.img(HMRFseg.nucleus_mask(:)), "FaceColor", [1, 1, 1]/6, "EdgeColor", "none");
% binwidth = h.BinWidth;
% for i = 1:nclust
% u = HMRFseg.mu(i);
% sig = HMRFseg.sigma(i);
% x = u-3*sig:0.001:u+3*sig;
% y = HMRFseg.a_counts(i)*normpdf(x, u, sig)*binwidth;
% plot(x, y, 'r', "LineWidth",1);
% end
% hold off
% xlabel("Intensity");
% ylabel("Counts");
% print(fig1, [save_path, filesep, filename(1:(end-4)), '_intensity.pdf'], '-dpdf');


% ----- Calculate enrichment of D in each chromatin class ------ %
logD_range = linspace(-3,1,10);
D_range = 10.^logD_range;

% Append additional edge values to classify all D < 0.001 and D > 10 correctly
D_range_edges = [0, D_range, Inf];

% Initialize counters for D range occurrence in each label
countDRangesInLabels = zeros(length(D_range_edges) - 1, nclust);

% Iterate over data points to classify them
for i = 1:numel(binX)
    if binX(i) > 0 && binY(i) > 0 && full_partition(binX(i), binY(i)) > 0
        % Get the label from full_partition
        label = full_partition(binX(i), binY(i));
        
        % Find which D range the current D value falls into
        dValue = tab_SR_pred.D(i);
        dRangeIndex = find(dValue >= D_range_edges(1:end-1) & dValue < D_range_edges(2:end), 1);
        
        % Increment the count for this D range and label
        if ~isempty(dRangeIndex)
            countDRangesInLabels(dRangeIndex, label) = countDRangesInLabels(dRangeIndex, label) + 1;
        end
    end
end

% Calculate enrichment for each D range in each label
totalInLabels = sum(countDRangesInLabels, 1); % Total count in each label
enrichment = bsxfun(@rdivide, countDRangesInLabels, totalInLabels);

% Display the results
fprintf('Enrichment of D ranges in each label:\n');
for dIdx = 1:size(enrichment, 1)
    fprintf('D Range [%.4f, %.4f):\n', D_range_edges(dIdx), D_range_edges(dIdx+1));
    for label = 1:numLabels
        fprintf('  Label %d: Count = %d, Enrichment = %.2f\n', ...
                label, countDRangesInLabels(dIdx, label), enrichment(dIdx, label) * 100);
    end
end

figure;
heatmap(enrichment,'Colormap',J2)



% -------------------------------------------------------------------
% Plot 3D line histograms for each chromatin class on vertical planes
% -------------------------------------------------------------------
figure;
hold on;

% Compute more finely spaced log-scaled bin edges
sharedEdges = linspace(-3, 2, 50); % Increased number of bins for finer resolution

colors = jet(nclust); % Generate distinct colors for each class

for k = 1:nclust
    if ~isempty(D_values{k})
        % Compute histogram counts with shared log-scaled bin edges
        counts = histcounts(log10(D_values{k}), sharedEdges);
        
        % Normalize counts to probabilities
        probabilities = counts / sum(counts);
        
        % Calculate bin centers in log space
        binCenters = (sharedEdges(1:end-1) + sharedEdges(2:end)) / 2;

        % Apply moving average to smooth probabilities
        windowSize = 5; % Define the window size for the moving average
        smoothProbabilities = movmean(probabilities, windowSize);

        % Plot as a line in 3D: X -> class index (k), Y -> log-scaled binCenters, Z -> probabilities
        plot3(binCenters, repmat(k, size(binCenters)), smoothProbabilities, ...
              'LineWidth', 2, 'Color', colors(k, :));
        
        % Create vertices for the polygon (including base at z=0)
        xVertices = [binCenters(1), binCenters, binCenters(end)];
        yVertices = [k,repmat(k, size(binCenters)), k];
        zVertices = [0,smoothProbabilities, 0];
        
        % Fill under the line plot
        fill3(xVertices, yVertices, zVertices, colors(k, :), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    end
end

% Labeling and view settings
xlabel('Log-scaled D Value');
ylabel('Chromatin Class ');
zlabel('Probability');
ax = gca; 
ax.YLim = [1, nclust];
%title('3D Line Histograms of D for Each Chromatin Class on Vertical Planes (Log Scale)');
view(3);
grid on;


% -------------------------------------------------------------------
% Plot 2D line histograms for each chromatin class on vertical planes
% -------------------------------------------------------------------
figure;
hold on;

% Compute more finely spaced log-scaled bin edges
sharedEdges = linspace(-3, 2, 50); % Increased number of bins for finer resolution

colors = jet(nclust); % Generate distinct colors for each class

for k = 1:nclust
    if ~isempty(D_values{k})
        % Compute histogram counts with shared log-scaled bin edges
        counts = histcounts(log10(D_values{k}), sharedEdges);
        
        % Normalize counts to probabilities
        probabilities = counts / sum(counts);
        
        % Calculate bin centers in log space
        binCenters = (sharedEdges(1:end-1) + sharedEdges(2:end)) / 2;

        % Apply moving average to smooth probabilities
        windowSize = 5; % Define the window size for the moving average
        smoothProbabilities = movmean(probabilities, windowSize);

        % Plot as a line in 3D: X -> class index (k), Y -> log-scaled binCenters, Z -> probabilities
        plot(binCenters, smoothProbabilities, ...
              'LineWidth', 2, 'Color', colors(k, :));
        
                
    end
end

% Labeling and view settings
xlabel('Log-scaled D Value');
ylabel('Probability ');
% Create a cell array of legend entries for chromatin classes 1 to 7
legendEntries = arrayfun(@(i) sprintf('Chromatin class %d', i), 1:7, 'UniformOutput', false);
% Add the legend to the plot
legend(legendEntries);
title('3D Line Histograms of D for Each Chromatin Class on Vertical Planes (Log Scale)');

% grid on;

%% Generate CD/non-IC boundary mask image
% Load or define your labeled image
% full_nucleus_partition = imread('your_image.png'); % Example if using an image file

% Assuming `full_nucleus_partition` is already defined and contains labels from 0 to 7

% Step 1: Scale the original image 5 times larger
scaleFactor = 2;
resizedImage = imresize(full_nucleus_partition, scaleFactor, 'nearest');
resize_CD = imresize(simple_nucleus_partition,scaleFactor,"nearest");
resize_immobileImPALM = imresize(immobileImPALM,scaleFactor,"nearest");


% Use edge detection on the upscaled region to get a boundary of one-pixel width
boundary_CD = bwperim(resize_CD);
resizedImage_nonIC = resizedImage > 1;
boundary_nonIC = bwperim(resizedImage_nonIC);


% save_path = 'C:\Users\Zuhui\Downloads\MPALM_paper_result\Figure4567_render_example\20231013_Clust02_EF1a-H2B_642_25uMPA646_30p5ms_20kframe_02_Cell01\UNet_mask_MBX_20240620_2035_epoch20_Ch1_SR_pred_v3_driftCorr\HMRF';
boundary_CD(boundary_CD>0) = 10000;
% imwrite(uint16(boundary_CD), [save_path, filesep, 'CD_boundary_5xScale.tif'], 'WriteMode', 'overwrite');
boundary_nonIC(boundary_nonIC>0) = 10000;
% imwrite(uint16(boundary_nonIC), [save_path, filesep, 'nonIC_boundary_5xScale.tif'], 'WriteMode', 'overwrite');
% imwrite(uint16(resize_immobileImPALM), [save_path, filesep, 'immobilePALM_5xScale.tif'], 'WriteMode', 'overwrite');


%% Plot PALM/Diffusion Map/MPALM for different D range
close all;
D_range = [0 0.05 0.5 1.5 Inf];

for i = 4% 1:4
    tab_SR_pred_filtD = tab_SR_pred(tab_SR_pred.D>=D_range(i) & tab_SR_pred.D<D_range(i+1),:);
    params.zoomROI = [567,281,104,107]; % x y w h
    
    params.dxy = 30; % rendered pixel size
    params.sigma_render = 50;
    params.camPixelSize_Ch1 = 110;
    
    params.Ch1_loc_XYLim = [0,params.camPixelSize_Ch1*params.ImWidth,0,params.camPixelSize_Ch1*params.ImHeight];
    
    params.render_minlogD = -2; % Min of log10(D)
    params.render_maxlogD = 0.699; % Max of log10(D)

    params.ColorscaleSwitch = "Log";
    switch params.ColorscaleSwitch
        case "Log"
            % log scale D colormap
            MinLogD = params.render_minlogD; % Min of log10(D)
            MaxLogD = params.render_maxlogD; % Max of log10(D)                    
            params.D_edge = 10.^linspace(MinLogD, MaxLogD, 256);

            % % Find the colormap mapping index of each range of logD. and
            % % map the first black color to all NaN data later
            % logDValues = linspace(MinLogD, MaxLogD, 256);
            % DValues = 10.^logDValues;
            % app.areaEdge = app.D_SRArea_const.Value*sqrt(DValues);
        case "Linear"
            % Linear scale D colormap
            MinLogD = params.render_minlogD; % Min of log10(D)
            MaxLogD = params.render_maxlogD; % Max of log10(D)         
            params.D_edge = linspace(10^MinLogD, 10^MaxLogD, 256);
            % MinD = 10^(app.render_minlogD.Value);
            % MaxD = 10^(app.render_maxlogD.Value);
            % DValues = linspace(MinD, MaxD, 256);
            % app.areaEdge = app.D_SRArea_const.Value*sqrt(DValues);
    end
   
    rgbColors = flipud(hsv(256));rgbColors = rgbColors(39:end,:);
    params.cmap_mpalm = interp1(linspace(0, 1, size(rgbColors, 1)), rgbColors, linspace(0, 1, 256));  % keep one bin index for NaN value;
    
    params.render_minDensity_Ch1 = 0;
    params.render_maxDensity_Ch1 = 5000;
    
    params.mpalm_DensityMin = 0;
    params.mpalm_DensityMax = 2000;

    
    
    
    [Ch1Img_filt_densityPALM,rgbImOut_mobility,img_mpalm] = MPALM_RENDER_forScript(tab_SR_pred_filtD,params);
    
    
    fig1 = figure('Position',[63 434 560 420],'Name',sprintf('D range between %g ~ %g',D_range(i),D_range(i+1)));
    ax_img_palm = axes(fig1,'Position',[0 0 1 1]);
    ax_CLim_PALM = [params.render_minDensity_Ch1,params.render_maxDensity_Ch1];
    imshow(Ch1Img_filt_densityPALM,ax_CLim_PALM,'Parent',ax_img_palm,Colormap=hot)
    ax_img_palm.XTick = [];ax_img_palm.YTick = [];

    % Set the limits to zoom into the ROI
    ax_img_palm.XLim = [params.zoomROI(1), params.zoomROI(1) + params.zoomROI(3)];
    ax_img_palm.YLim = [params.zoomROI(2), params.zoomROI(2) + params.zoomROI(4)];        
    
    fig2 = figure('Position',[627 434 560 420],'Name',sprintf('D range between %g ~ %g',D_range(i),D_range(i+1)));
    ax_img_mobility = axes(fig2,'Position',[0 0 1 1]);
    imshow(rgbImOut_mobility,'Parent',ax_img_mobility)
    ax_img_mobility.XTick = [];ax_img_mobility.YTick = [];
    % Set the limits to zoom into the ROI
    ax_img_mobility.XLim = [params.zoomROI(1), params.zoomROI(1) + params.zoomROI(3)];
    ax_img_mobility.YLim = [params.zoomROI(2), params.zoomROI(2) + params.zoomROI(4)];
    
    fig3 = figure('Position',[1191 435 560 420],'Name',sprintf('D range between %g ~ %g',D_range(i),D_range(i+1)));
    ax_CLim_MPALM = axes(fig3,'Position',[0 0 1 1]);
    imshow(img_mpalm,'Parent',ax_CLim_MPALM);
    ax_CLim_MPALM.XTick = [];ax_CLim_MPALM.YTick = [];
    % Set the limits to zoom into the ROI
    ax_CLim_MPALM.XLim = [params.zoomROI(1), params.zoomROI(1) + params.zoomROI(3)];
    ax_CLim_MPALM.YLim = [params.zoomROI(2), params.zoomROI(2) + params.zoomROI(4)];    
    
    % ---- Overlay CD, non-IC boundary ---- %
    % Normalize the image data according to the specified axes limits
    normalizedImage = mat2gray(Ch1Img_filt_densityPALM, ax_CLim_PALM);   
    % Apply the colormap 'hot' to the normalized image
    colormapData = colormap(gray(256)); % Generate a 256-level hot colormap
    rgbImage = ind2rgb(im2uint8(normalizedImage), colormapData);
    % Convert the RGB image to uint8
    rgbImage = im2uint8(rgbImage);
    PlotBoundaryOverlayImg(rgbImage,boundary_CD,boundary_nonIC,scaleFactor,params.zoomROI);
    
    rgbImage = im2uint8(rgbImOut_mobility);
    PlotBoundaryOverlayImg(rgbImage,boundary_CD,boundary_nonIC,scaleFactor,params.zoomROI);

    rgbImage = im2uint8(img_mpalm);
    PlotBoundaryOverlayImg(rgbImage,boundary_CD,boundary_nonIC,scaleFactor,params.zoomROI);
end



%% Plot PALM/MPALM/Diffusion map/HMRF segmentation results in the zoomed ROI using all locs
close all;
PALM_img = imread('C:\Users\Zuhui\Downloads\MPALM_paper_result\Figure4567_render_example\20231013_Clust02_EF1a-H2B_642_25uMPA646_30p5ms_20kframe_02_Cell01\UNet_mask_MBX_20240620_2035_epoch20_Ch1_SR_pred_v3_driftCorr\PALM_Ch1_gray.tif');
MPALM_img = imread('C:\Users\Zuhui\Downloads\MPALM_paper_result\Figure4567_render_example\20231013_Clust02_EF1a-H2B_642_25uMPA646_30p5ms_20kframe_02_Cell01\UNet_mask_MBX_20240620_2035_epoch20_Ch1_SR_pred_v3_driftCorr\MPALM_new.tif');
DiffMap_img = imread('C:\Users\Zuhui\Downloads\MPALM_paper_result\Figure4567_render_example\20231013_Clust02_EF1a-H2B_642_25uMPA646_30p5ms_20kframe_02_Cell01\UNet_mask_MBX_20240620_2035_epoch20_Ch1_SR_pred_v3_driftCorr\MobilityMap_new.tif');

PlotBoundaryOverlayImg(PALM_img,boundary_CD,boundary_nonIC,scaleFactor,params.zoomROI);
PlotBoundaryOverlayImg(MPALM_img,boundary_CD,boundary_nonIC,scaleFactor,params.zoomROI);
PlotBoundaryOverlayImg(DiffMap_img,boundary_CD,boundary_nonIC,scaleFactor,params.zoomROI);

% Plot HMRF results in zoomed ROI
SE = strel("disk",7);
temp_full_nucleus_partition = full_nucleus_partition;
erodedBW = imerode(transpose(HMRFseg.nucleus_mask),SE);
temp_full_nucleus_partition(~erodedBW) = 0;

figure; ax = axes('Position',[0 0 1 1]);
imagesc(ax,temp_full_nucleus_partition)
colormap(ax,cmap_all)
axis(ax,'off')
axis image ij
% Set the limits to zoom into the ROI
ax.XLim = [params.zoomROI(1), params.zoomROI(1) + params.zoomROI(3)];
ax.YLim = [params.zoomROI(2), params.zoomROI(2) + params.zoomROI(4)]; 


% %% Wangbo's IC/PC/CD boundary plot
% 
% save_path = 'C:\Users\Zuhui\Downloads\MPALM_paper_result\Figure4567_render_example\20231013_Clust02_EF1a-H2B_642_25uMPA646_30p5ms_20kframe_02_Cell01\UNet_mask_MBX_20240620_2035_epoch20_Ch1_SR_pred_v3_driftCorr\HMRF';
% 
% % find local maximum to define chromatin domain center
% tempPALM = reshape(HMRFseg.img, size(HMRFseg.nucleus_mask)); ImPALM = tempPALM';
% CDmask = reshape(HMRFseg.img_class > 3, size(HMRFseg.nucleus_mask));
% tempPALM = ImPALM; tempPALM(~HMRFseg.nucleus_mask') = 0;
% tempPALM(~CDmask') = 0;
% CDcenter = imregionalmax(tempPALM, 8);
% CDcenter = uint16(CDcenter);
% CDcenter(CDcenter>0) = 10000;
% tempCDcenter = imdilate(CDcenter, strel('disk', 1));
% imwrite(uint16(tempCDcenter), [save_path, filesep, 'CDcenter.tif'], 'WriteMode', 'overwrite');
% 
% % find IC interface:defined as pixel of class2 with 8-connected
% % nearby class1 pixel
% temp_img_class = zeros(size(HMRFseg.nucleus_mask, 1)+2, size(HMRFseg.nucleus_mask, 2)+2);
% temp_img_class(2:(end-1), 2:(end-1)) = reshape(HMRFseg.img_class, size(HMRFseg.nucleus_mask));
% ICinterface = zeros(size(HMRFseg.nucleus_mask));
% for x = 2:(size(temp_img_class, 2)-1)
%     for y = 2:(size(temp_img_class, 1)-1)
%         if temp_img_class(y, x) == 2
%             temp_neighbor = temp_img_class((y-1):(y+1), (x-1):(x+1));
%             if sum(temp_neighbor(:) == 1) > 0
%             ICinterface(y-1, x-1) = 10000;
%             end
%         end
%     end
% end
% ICinterface = ICinterface';
% imwrite(uint16(ICinterface), [save_path, filesep, filename1(1:(end-4)), '_ICinterface.tif'], 'WriteMode', 'overwrite');
% 
% % find CD interface:defined as pixel of class4 with 8-connected
% % nearby class3 pixel
% temp_img_class = zeros(size(HMRFseg.nucleus_mask, 1)+2, size(HMRFseg.nucleus_mask, 2)+2);
% temp_img_class(2:(end-1), 2:(end-1)) = reshape(HMRFseg.img_class, size(HMRFseg.nucleus_mask));
% CDinterface = zeros(size(HMRFseg.nucleus_mask));
% for x = 2:(size(temp_img_class, 2)-1)
%     for y = 2:(size(temp_img_class, 1)-1)
%         if temp_img_class(y, x) == 3+1
%             temp_neighbor = temp_img_class((y-1):(y+1), (x-1):(x+1));
%             if sum(temp_neighbor(:) == 3) > 0
%                 CDinterface(y-1, x-1) = 10000;
%             end
%         end
%     end
% end
% CDinterface = CDinterface';
% imwrite(uint16(CDinterface), [save_path, filesep, filename1(1:(end-4)), '_CDinterface.tif'], 'WriteMode', 'overwrite');

%% Auxiliary Function
function PlotBoundaryOverlayImg(img,boundary_CD,boundary_nonIC,scaleFactor,varargin)
%SaveBoundaryOverlayImg Summary of this function goes here
%   Detailed explanation goes here

if size(varargin) == 1
    has_roi = true;
    save_img = false;
    roi_coord = varargin{1}*scaleFactor;
elseif size(varargin) == 2
    has_roi = false;
    save_img = true;
    save_path = varargin{1};
    save_name = varargin{2};
elseif size(varargin) == 3
    has_roi = true;
    save_img = true;
    roi_coord = varargin{1}*scaleFactor;
    save_path = varargin{2};
    save_name = varargin{3};
end


resized_img = imresize(img, scaleFactor, 'nearest');

% Ensure both boundary masks are logical
boundary_CD = logical(boundary_CD);
boundary_nonIC = logical(boundary_nonIC);

% Create a copy of the resized RGB image for overlay
overlayImage = resized_img;

% Set boundary_CD area to yellow [R, G, B] = [255, 255, 0]
overlayImage(boundary_CD) = 255; % Red channel
overlayImage(:, :, 2) = max(overlayImage(:, :, 2), uint8(boundary_CD) * 255); % Green channel
overlayImage(:, :, 3) = min(overlayImage(:, :, 3), uint8(~boundary_CD) * 255); % Blue channel remains 0 where mask is true

% Now set boundary_nonIC area to white [R, G, B] = [255, 255, 255]
overlayImage(boundary_nonIC) = 255; % Red channel
overlayImage(:, :, 2) = max(overlayImage(:, :, 2), uint8(boundary_nonIC) * 255); % Green channel
overlayImage(:, :, 3) = max(overlayImage(:, :, 3), uint8(boundary_nonIC) * 255); % Blue channel

% Display the overlaid image with both boundaries
figure('Visible','on'); ax = axes('Position',[0 0 1 1]);
imshow(overlayImage,'Parent',ax);
if has_roi
    ax.XTick = [];ax.YTick = [];
    % Set the limits to zoom into the ROI
    ax.XLim = [roi_coord(1), roi_coord(1) + roi_coord(3)];
    ax.YLim = [roi_coord(2), roi_coord(2) + roi_coord(4)];  
end
if save_img
    imgObj = get(ax, 'Children');            
    % Get the image data from the CData property
    imageData = imgObj.CData;
    % Save the image data using imwrite
    imwrite(imageData, fullfile(save_path, sprintf('%s.tif',save_name)));
end
end
