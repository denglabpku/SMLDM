function [Xgrid,Ygrid,Zgrid] = ARRAY2GRID(varargin)
%ARRAY2GRID Convert PSF array into grid data point.
%   If 1 input: xy, use raw Zgrid (intensity) to draw raw PSF
%   profile.
%   If 2 input: xyI and BoundxyI, fill the NaN Zgrid by noise from boundary pixel 

%% Convert individual PSF xyI stored in cell_PSF into matrix
xyI = varargin{1};
% xyI = cell_PSF.xyI{1}{1};
X = xyI(:,1); Y = xyI(:,2); Z = xyI(:,3);
Xs = unique(X); % image X coordinate of pixel of PSF
Ys = unique(Y); % image Y coordinate of pixel of PSF
Xi = arrayfun(@(x) find(Xs == x), X); %return index of Xs that equal to X
Yi = arrayfun(@(y) find(Ys == y), Y); %return index of Ys that equal to Y
Li = Yi + (Xi -1) * numel(Ys);
Zgrid = nan(numel(Ys), numel(Xs));
Zgrid(Li) = Z; %pixel intensity with same shape as original image

[Xgrid,Ygrid] = meshgrid(Xs,Ys);

if nargin == 2
    %% Substract background noise
    BoundxyI = varargin{2};
    
    % Int_noise = mean(cell_PSF.BoundxyI{i_plane}{PSF_idx}(:,3)); % calculate intensity noise as the mean of PSF_boundLabel pixel intensity
    % Int_noise = mean(BoundxyI(:,3));

    % BoundxyI = cell_PSF.BoundxyI{1}{1}
    Int_noise = mean(BoundxyI(:,3));
    Zgrid(isnan(Zgrid)) = Int_noise; % Fill the NaN of PSF rectangular matrix with noise
    % Zgrid = Zgrid-Int_noise; % Substract the noise from each pixel intensity
    % Zgrid(Zgrid < 0) = 0; % make the negative intensity to zero
end
end

