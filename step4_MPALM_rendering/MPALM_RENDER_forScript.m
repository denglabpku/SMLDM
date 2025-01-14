function [Ch1Img_filt_densityPALM,rgbImOut_mobility,img_mpalm] = MPALM_RENDER_forScript(loc_table,params)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% ================== 2023/09/21 mobility-PALM rendering without cluster analysis ===================== %%
%% MPALM render considering distance and local density followed by cross-one non-NaN convolution
dxy = params.dxy; % rendered pixel size, unit NM
sigma_render = params.sigma_render;
Xmin = params.Ch1_loc_XYLim(1);Xmax=params.Ch1_loc_XYLim(2);Ymin = params.Ch1_loc_XYLim(3);Ymax=params.Ch1_loc_XYLim(4);
Edges{1}=Xmin:dxy:Xmax+dxy; % make edge cover all data points
Edges{2}=Ymin:dxy:Ymax+dxy;
             
%% Step2: Generate mobility map with proper gaussian weighted average
Xpos = loc_table.Xpos;
Ypos = loc_table.Ypos;
D = loc_table.D;
        
[Im,Xedges,Yedges,binX,binY] = histcounts2(Xpos,Ypos,Edges{1},Edges{2}); % bin every localizations into grid, similar to Im = hist3([Xpos',Ypos'],'Edges',Edges);

% remove points not fit the grid
Xpos = Xpos(~binX==0);
Ypos = Ypos(~binY==0);
D = D(~binX==0);
[Im,Xedges,Yedges,binX,binY] = histcounts2(Xpos,Ypos,Edges{1},Edges{2});

% Calculate the average Cdata for each bin index
matrix = accumarray([binX(:), binY(:)], D(:), size(Im), @mean);

% Im = hist3([Xpos',Ypos'],'Edges',Edges);

TempX=-round(3*sigma_render/dxy)*dxy:dxy:round(3*sigma_render/dxy)*dxy;
TempY=-round(3*sigma_render/dxy)*dxy:dxy:round(3*sigma_render/dxy)*dxy;

ConVecX = exp(-0.5*(TempX/sigma_render).^2);
ConVecX=ConVecX/sum(ConVecX);
ConVecY = exp(-0.5*(TempY/sigma_render).^2);
ConVecY=ConVecY/sum(ConVecY);

Im_density = conv2(ConVecX,ConVecY,Im,'same');

Im2=Im_density/dxy/dxy*10^6; % convert unit of localizations in each bin from nm^-2 to um^-2
Ch1Img_filt_densityPALM = transpose(Im2); % x-y direction of histcounts2 are transposed with original image

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
% ImOut_mobility = nanconv(filteredMatrix, kernel,"nanout"); % this get sharper results
ImOut_mobility = nanconv(maskedA, kernel,"nanout"); % this get sharper results            
% smoothedA = nanconv(maskedA, kernel);

img_mobility = transpose(ImOut_mobility);
fprintf('mPALM render done \n')

%% Step3: combine PALM and mobility map

% log scale D colormap
MinLogD = params.render_minlogD; % Min of log10(D)
MaxLogD = params.render_maxlogD; % Max of log10(D)

rescale_ImOut_mobility = img_mobility;

% areaEdge = app.areaEdge;
rescale_ImOut_mobility(rescale_ImOut_mobility < 10^MinLogD & rescale_ImOut_mobility > 0 ) = 10^MinLogD;
rescale_ImOut_mobility(rescale_ImOut_mobility >  10^MaxLogD) = 10^MaxLogD;

% Create custom colormap, first black reserved for NaN
J = params.cmap_mpalm;
J(1,:) = [0 0 0];
local_cmap_mpalm = J;

CData_rescale_ImOut_mobility = discretize(rescale_ImOut_mobility,params.D_edge);             
CData_rescale_ImOut_mobility = CData_rescale_ImOut_mobility +1;
CData_rescale_ImOut_mobility(isnan(CData_rescale_ImOut_mobility)) = 1;
%             hI.CData = CData_rescale_ImOut_mobility;
%             hI.CDataMapping = 'direct';                                   
%             ax = gca; ax.XTickLabel = '';ax.YTickLabel = '';
%             ax.Units = 'normalized';ax.Position = [0 0 1 1];

%---------- PLOT: PALM density integrated MPALM using HSV color scheme ------%
rgbImOut_mobility = ind2rgb(CData_rescale_ImOut_mobility,J);




          

hsv_image = rgb2hsv(rgbImOut_mobility);
% lets construct the HSV coded mobility-PALM image
ax_CLim = [params.mpalm_DensityMin params.mpalm_DensityMax];
rescale_ax_CLim = [ax_CLim(1)/max(Ch1Img_filt_densityPALM,[],'all'), ax_CLim(2)/max(Ch1Img_filt_densityPALM,[],'all')];
if min(rescale_ax_CLim) < 0 || max(rescale_ax_CLim) > 1
    rescale_ax_CLim = [0 1];
end
rescale_ImOut_density = rescale(Ch1Img_filt_densityPALM);
scaled_ImOut_density = imadjust(rescale_ImOut_density,rescale_ax_CLim);

% >>>> Option1: black background for MPALM
hsv_image(:,:,3) = scaled_ImOut_density; % Value for PALM localization density
% <<<<

% >>>> Option2: white background for MPALM
% app.hsv_image(:,:,2) = scaled_ImOut_density;
% app.hsv_image(:,:,3) = ones(size(scaled_ImOut_density));
% <<<<

img_mpalm = hsv2rgb(hsv_image);

end