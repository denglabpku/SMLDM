function [ImOut_density, img_mpalm,img_density_mpalm, fig1,fig2,fig3] = MPALM_RENDERv3(raw_table,RenderParam)

%MPALM_RENDERv2 Summary of this function goes here
%   This function will render the density-PALM and mobility map and merge them together.
%   When merging, the colormap uses HSV color space, where the hue represent the mobility,
%   the saturation indicates the pixel contains localization(signal) or not, and the value represents
%   the density (regular PALM). When mapping the mobility (area of motion-blur) to hue, the mobility map
%   will first scaled so that only area >= 50 will be rendered, and the area from 40 to 80 will be linearly
%   mapped to blue (hue 0) and red (hue 0.66), area higher than 80 will be mapped to red (hue 0.66).

% % Debug purpose only
% raw_table = NPWAD3_tab_psf_fitresult(NPWAD3_tab_psf_fitresult.TotalPhoton <= 2000,:); % filter out abnormal molecules
% raw_table.Xpos = raw_table.Xpos * 110;
% raw_table.Ypos = raw_table.Ypos * 110;

% RenderParam.dxy = 50;
% RenderParam.sigma_render = 50;
% RenderParam.AreaCutoff = [50 80];
% RenderParam.DensityCutoff = [0 500];
% RenderParam.Xrange = [min(raw_table.Xpos) max(raw_table.Xpos)];
% RenderParam.Yrange = [min(raw_table.Ypos) max(raw_table.Ypos)];

%% ================== 2023/08/17 newer version of mobility-PALM rendering without cluster analysis ===================== %%
D_SRArea_const = RenderParam.D_SRArea_const;
flipped = true;
dxy = RenderParam.dxy; % rendered pixel size 20 nm
sigma_render = RenderParam.sigma_render;
Xmin = RenderParam.Xrange(1);Xmax=RenderParam.Xrange(2);Ymin = RenderParam.Yrange(1);Ymax=RenderParam.Yrange(2);
Edges{1}=Xmin:dxy:Xmax+dxy; % make edge cover all data points
Edges{2}=Ymin:dxy:Ymax+dxy;

Xpos = raw_table.Xpos; % get immobile molecules info
Ypos = raw_table.Ypos;
Area = raw_table.SRArea;

%% Step1: generate density based PALM image using immobile molecules
[Im,Xedges,Yedges,binX,binY] = histcounts2(Xpos,Ypos,Edges{1},Edges{2}); % bin every localizations into grid, similar to Im = hist3([Xpos',Ypos'],'Edges',Edges);

TempX=-round(3*sigma_render/dxy)*dxy:dxy:round(3*sigma_render/dxy)*dxy;
TempY=-round(3*sigma_render/dxy)*dxy:dxy:round(3*sigma_render/dxy)*dxy;

ConVecX = exp(-0.5*(TempX/sigma_render).^2);
ConVecX=ConVecX/sum(ConVecX);
ConVecY = exp(-0.5*(TempY/sigma_render).^2);
ConVecY=ConVecY/sum(ConVecY);
% ZW's updated shorten version, which gives same results
Im2 = conv2(ConVecX,ConVecY,Im,'same');
Im2=Im2/dxy/dxy*10^6; % convert unit of localizations in each bin from nm^-2 to um^-2
ImOut_density = transpose(Im2); % x-y direction of histcounts2 are transposed with original image

%% Step2: Generate mobility map with proper gaussian weighted average
% Calculate the average Cdata for each bin index
matrix = accumarray([binX(:), binY(:)], Area(:), size(Im), @mean);

TempX=-round(3*sigma_render/dxy)*dxy:dxy:round(3*sigma_render/dxy)*dxy;
TempY=-round(3*sigma_render/dxy)*dxy:dxy:round(3*sigma_render/dxy)*dxy;

ConVecX = exp(-0.5*(TempX/sigma_render).^2);
ConVecX=ConVecX/sum(ConVecX);
ConVecY = exp(-0.5*(TempY/sigma_render).^2);
ConVecY=ConVecY/sum(ConVecY);

Im_density = conv2(ConVecX,ConVecY,Im,'same');

% Create a logical mask indicating the zero elements
zeroMask = (matrix == 0);

% Compute the weighted average of the Gaussian distribution for zero elements
if any(zeroMask(:))

    % Get the size of the matrix
    [M, N] = size(matrix);
    
    % Initialize the filtered matrix
    filteredMatrix = matrix;
    
    % Loop over each element in the matrix
    for i = 1:M
        for j = 1:N
            if zeroMask(i, j) % if this pixel is zero
                % Initialize total weight and weighted sum
                totalWeight = 0;
                weightedSum = 0;
                
                % Loop over neighbors (a 3x3 pixel square center by the current zero-value pixel when sigma_render equals dxy)
                for dx = -round(1*sigma_render/dxy):round(1*sigma_render/dxy)  %-round(2*sigma_render/dxy):round(2*sigma_render/dxy)
                    for dy = -round(1*sigma_render/dxy):round(1*sigma_render/dxy) %-round(2*sigma_render/dxy):round(2*sigma_render/dxy)
                        % Check if neighbor is within the matrix bounds
                        if i+dx >= 1 && i+dx <= M && j+dy >= 1 && j+dy <= N
                            % Check if neighbor is non-zero
                            if ~zeroMask(i+dx, j+dy)
                                % Compute weights based on Gaussian distribution
                                weightX = exp(-0.5*(dx*dxy/sigma_render).^2); %exp(-0.5*(dx)^2);
                                weightY = exp(-0.5*(dy*dxy/sigma_render).^2); %exp(-0.5*(dy)^2);
                                
                                % Accumulate weighted sum and total weight
                                weightedSum = weightedSum + matrix(i+dx, j+dy) * weightX * weightY * Im_density(i+dx, j+dy);
                                totalWeight = totalWeight + weightX * weightY * Im_density(i+dx, j+dy);
                            end
                        end
                    end
                end
                
                % Compute the weighted average for the zero element
                filteredMatrix(i, j) = weightedSum / totalWeight;
            end
        end
    end
else
    % If there are no zero elements, the filtered matrix is the same as the original matrix
    filteredMatrix = matrix;
end

% Create a binary mask to identify the non-zero values in A
mask = filteredMatrix ~= 0;

% % Define the kernel mask for smoothing, not as good as
% cross-on kernel
% kernel = [1, 1, 1;
%           1, 1, 1;
%           1, 1, 1];

% cross-one kernel
kernel = [0, 1, 0;
          1, 1, 1;
          0, 1, 0];

% Normalize the kernel so that its sum is equal to 1
kernel = kernel / sum(kernel(:));

% Create a matrix with NaN values at zero positions and non-zero values at non-zero positions in A
maskedA = filteredMatrix;
maskedA(mask == 0) = NaN;

% Perform kernel smoothing on the non-zero elements in A using nanconv
ImOut_mobility = nanconv(maskedA, kernel,"nanout"); % this get sharper results
% smoothedA = nanconv(maskedA, kernel);

img_mpalm = transpose(ImOut_mobility);

%% Step3: combine PALM and mobility
% ==============  set colormap for mPALM render ============== %
rgbColors = flipud(hsv(256));
rgbColors = rgbColors(86:end,:);
cmap_mpalm = interp1(linspace(0, 1, size(rgbColors, 1)), rgbColors, linspace(0, 1, 256));

% log scale D colormap
MinLogD = RenderParam.render_minlogD; % Min of log10(D)
MaxLogD = RenderParam.render_maxlogD; % Max of log10(D)

% Find the colormap mapping index of each range of logD. and
% map the first black color to all NaN data later
logDValues = linspace(MinLogD, MaxLogD, 256);
DValues = 10.^logDValues;
areaEdge = D_SRArea_const*sqrt(DValues);

%------------PLOT: MPALM and its copied RGB image ------------%
rescale_ImOut_mobility = img_mpalm;
rescale_ImOut_mobility(rescale_ImOut_mobility < min(areaEdge) & rescale_ImOut_mobility > 0 ) = min(areaEdge);
rescale_ImOut_mobility(rescale_ImOut_mobility >  max(areaEdge)) = max(areaEdge);

% Create custom colormap, first black reserved for NaN
J = cmap_mpalm;
J(1,:) = [0 0 0];

CData_rescale_ImOut_mobility = discretize(rescale_ImOut_mobility,areaEdge); 
CData_rescale_ImOut_mobility = CData_rescale_ImOut_mobility +1;
CData_rescale_ImOut_mobility(isnan(CData_rescale_ImOut_mobility)) = 1;

%---------- PLOT: PALM density integrated MPALM using HSV color scheme ------%
rgbImOut_mobility = ind2rgb(CData_rescale_ImOut_mobility,J);

hsv_image = rgb2hsv(rgbImOut_mobility);
% lets construct the HSV coded mobility-PALM image
ax_CLim = RenderParam.DensityCutoff_densityMPALM;
rescale_ImOut_density = rescale(ImOut_density);
scaled_ImOut_density = imadjust(rescale_ImOut_density,[ax_CLim(1)/max(ImOut_density,[],'all'), ax_CLim(2)/max(ImOut_density,[],'all')]);

% >>>> Option1: black background for MPALM
hsv_image(:,:,3) = scaled_ImOut_density; % Value for PALM localization density
% <<<<

% >>>> Option2: white background for MPALM
% app.hsv_image(:,:,2) = scaled_ImOut_density;
% app.hsv_image(:,:,3) = ones(size(scaled_ImOut_density));
% <<<<

img_density_mpalm = hsv2rgb(hsv_image);

%%
% scale ImOut_density linearly into [0, 1] using imcontrast is more easier
fig1 = figure('Visible','on','Units','inches','Position',[2.1 1.2 10 7]);
ax1 = axes(...
    'Units','normalized',...
    'Position', [0 0 1 1],... % stretch to fill axes 
    'XTickLabel','',...
    'YTickLabel','');
imshow(ImOut_density,RenderParam.DensityCutoff_densityPALM,'InitialMagnification',1000,'colormap',hot,'Parent',ax1);

fig2 = figure('Visible','on','Units','inches','Position',[2.1 1.2 10 7]);
ax2 = axes(...
    'Units','normalized',...
    'Position', [0 0 1 1],... % stretch to fill axes 
    'XTickLabel','',...
    'YTickLabel','');
imshow(rgbImOut_mobility,'InitialMagnification',1000,'Parent',ax2)

fig3 = figure('Visible','on','Units','inches','Position',[2.1 1.2 10 7]);
ax3 = axes(...
    'Units','normalized',...
    'Position', [0 0 1 1],... % stretch to fill axes 
    'XTickLabel','',...
    'YTickLabel','');
imshow(img_density_mpalm,'InitialMagnification',1000,'Parent',ax3)

end