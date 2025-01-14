function [fitresult, resnorm, exitflag] = PARA_GAUSS2DFIT(cell_PSF_xyI,cell_PSF_BoundxyI)
%PARA_GAUSS2DFIT Parallel fitting extracted PSF with general elliptical 2D
%Gaussian function
%   Detailed explanation goes here.
%   Dependency: Instrument Control Toolbox, parforprogress package to show
%   parallel computing progress.

%% Zuhui Wang
%% 2020/11/15
%% Convert individual PSF xyI stored in cell_PSF into matrix
[Xgrid,Ygrid,Zgrid] = cellfun(@ARRAY2GRID,cell_PSF_xyI,cell_PSF_BoundxyI,'UniformOutput',false);

% Fill the gap to match maxSz size
[Xgrid,Ygrid,Zgrid] = GRID_REGULARIZATION(Xgrid,Ygrid,Zgrid,cell_PSF_BoundxyI);


[xData, yData, zData] = cellfun(@prepareSurfaceData,Xgrid,Ygrid,Zgrid,'UniformOutput',false);


xData = xData.';
xData = cell2mat(xData);

xData = xData.';

yData = yData.';
yData = cell2mat(yData);
yData = yData.';

%% Set up the startpoint
% find starting points
[amp,ind] = cellfun(@max,zData);% amp is the amplitude,find the max Z as guess

numFit = size(Xgrid,1);
xo = NaN(numFit,1); %initial guess of x with max Z
yo = NaN(numFit,1); %initial guess of y with max Z
for i = 1:numFit
    xo(i) = xData(i,ind(i));
    yo(i) = yData(i,ind(i));
end

zData = zData.';
zData = cell2mat(zData);
zData = zData.';

ang = repmat(45,numFit,1); % angle in degrees.
sy = ones(numFit,1);
sx = ones(numFit,1);
zo = median(zData,2)-std(zData,0,2);
% zo = mean(zData,2);
xmax = max(xData,[],2)+2;
ymax = max(yData,[],2)+2;
xmin = min(xData,[],2)-2;
ymin = min(yData,[],2)-2;

% Set up fittype and options.(currently missed)

% Lower = [0, 0, 0, 0, xmin, ymin, 0];
Lower = zeros(numFit, 7);
Lower(:,[5,6]) = [xmin, ymin];

% Upper = [Inf, 180, Inf, Inf, xmax, ymax, Inf]; % angles greater than 90 are redundant
Upper = zeros(numFit,7);
Upper(:,[1,3,4,7]) = Inf(numFit,4);
Upper(:,2) = repmat(180,numFit,1);
Upper(:,[5,6]) = [xmax,ymax];

par = [amp, ang, sx, sy, xo, yo, zo];


xy(:,:,1) = xData;
xy(:,:,2) = yData;

tols = 1e-6; %tols = 1e-16;
options = optimoptions('lsqcurvefit',...
    'Algorithm','trust-region-reflective',...  changed from 'levenberg-marquardt' to 'trust-region-reflective'
    'Display','off',...
    'StepTolerance',tols,...
    'FunctionTolerance',tols);
%     'MaxFunctionEvaluations',7e2,...
%     'MaxIterations',7e2,...

fitresult = NaN(numFit,7);
resnorm = NaN(numFit,1);
exitflag = NaN(numFit,1);

% parfor_progress(numFit);
parfor i = 1:numFit
[fitresult(i,:),resnorm(i,:),~,exitflag(i,:),~] = ... fitresult(i,:) = ... %
    lsqcurvefit(@gaussian2D,par(i,:),xy(i,:,:),zData(i,:),Lower(i,:),Upper(i,:),options);
% parfor_progress;
end
% parfor_progress(0);

end
