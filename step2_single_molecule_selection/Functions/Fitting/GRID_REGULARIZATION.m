function [Xgrid,Ygrid,Zgrid] = GRID_REGULARIZATION(Xgrid,Ygrid,Zgrid,varargin)
%GRID_REGULARIZATION Regularize grid PSF data points to same dimension
%   Regularize grid PSF data points to same dimension. Append zero to
%   the output grid when fitting (varargin empty), or append NaN to the
%   output grid when plot fitting in FIT_3DIMG (varargin is true).

%% Zuhui Wang
%% 2020/11/15
%%
gridSz = cell2mat(cellfun(@size,Xgrid,'UniformOutput',false)); % [row#, col#]
maxSz = max(gridSz,[],1);%%%check
BoundxyI = varargin{1};

for i = 1: size(Xgrid,1)
if size(Xgrid{i},1) < maxSz(1) || size(Xgrid{i},2) < maxSz(2)
    
    maxX = maxSz(2) - size(Xgrid{i},2) + Xgrid{i}(1,end);
    newX = Xgrid{i}(1,1):1: maxX;
    
    maxY = maxSz(1) - size(Ygrid{i},1) + Ygrid{i}(end,1);
    newY = Ygrid{i}(1,1):1:maxY;
    
    [Xgrid{i},Ygrid{i}] = meshgrid(newX,newY);
    if isempty(varargin)
    % Extend Zgrid with 0, since the Zgrid is the noise subtracted, assume
    % append pixel are all no noise.
        temp1 = zeros(maxSz(1) - size(Zgrid{i},1),size(Zgrid{i},2));
        Zgrid{i} = vertcat(Zgrid{i},temp1);
        temp2 = zeros(size(Zgrid{i},1), maxSz(2) - size(Zgrid{i},2));
        Zgrid{i} = horzcat(Zgrid{i},temp2);
    else
%         % For plot Fit 3D only
%         temp1 = NaN(maxSz(1) - size(Zgrid{i},1),size(Zgrid{i},2));
%         Zgrid{i} = vertcat(Zgrid{i},temp1);
%         temp2 = NaN(size(Zgrid{i},1), maxSz(2) - size(Zgrid{i},2));
%         Zgrid{i} = horzcat(Zgrid{i},temp2);
        
        % Extend Zgrid with the mean of boundary pixel value
        Int_noise = mean(BoundxyI{i}(:,3));
        temp1 = Int_noise.*ones(maxSz(1) - size(Zgrid{i},1),size(Zgrid{i},2));
        Zgrid{i} = vertcat(Zgrid{i},temp1);
        temp2 = Int_noise.*ones(size(Zgrid{i},1), maxSz(2) - size(Zgrid{i},2));
        Zgrid{i} = horzcat(Zgrid{i},temp2);
    end
end
end

