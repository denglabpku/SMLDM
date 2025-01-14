function [corrs params] = Spatial_Temporal_Xcorr(data1, data2, mask, calib, dr,max_tau)
% corrs = Spatial_Temporal_Xcorr(data1, data2, mask, calib, dr,max_tau)
%
% SPATIAL_TEMPORAL_XCORR calculates pair-cross correlation between two sets of particles
% (data1 and data2) falling within a region of interest specified by mask, 
% and returns C(r,tau) in the structure CORRS, with fields C_11 and dC_11
% for the autocorrelation (G(r,tau)) of the DATA1, C_22 and dC_22 for the
% G(r,tau) of DATA2, and C_12 and dC_12 for the cross correlation of DATA1
% and DATA2. C_12 correlates DATA1(t) with DATA2(t+tau), whereas C_21
% correlates DATA2(t) with DATA1(t+tau). SPATIAL_TEMPORAL_XCORR also
% returns the structure params which contains the field radii, the
% interparticle distance for the bins C(r,tau) and G(r,tau) given in
% nanometer units. 
%
% DATA1 and DATA2 are structures with numel(DATA1) == numel(DATA2) = number 
% of frames to be correlated and have fields x and y corresponding to the x and y
% positions of particles in each frame. DATA1(i).x is a column of all x positions in
% frame i. DATA1(i).y is the column of y positions and is the same size as
% DATA1(i).x.  DATA1 and DATA2 should be in register with each other. 
%
% x and y should be in units of pixels, do not
% have to be integers, and should be given as column vectors
% for example [DATA1(i).x  DATA1(i).y] returns a 2 column matrix that gives
% the x and y position of all particles in frame i. 
%
% MASK is a logical matrix with 1's specifying the region of interest,
% having the same pixel units as DATA. size(MASK) = [xdim ydim]; thus the
% first dimension of MASK corresponds to DATA.x 
%
% CALIB is the length of one pixel in units of nanometers
%
% DR is the one dimensional bin length for the radii histogram M(r), in units of nanometers
%
% MAX_TAU is the maximum number of frames to correlate over. Set MAX_TAU to
% zero to only calculate C(r,tau=0). When MAX_TAU=0 G(r,tau) is not
% calculated for either channel. 
%
% NOTE: Autocorrelations can be normalized
% to probability density (PDF) by using Go(r,tau) = G(t,tau) - 1, and the
% PDF(r,tau) = Go(r,tau)/sum(Go(r,tau)*dr*mean_rads*2*pi)
%
%
% This code is part of the Supplemental Information for the revised manuscript: 
% Steady-state cross-correlations for live two-color super-resolution localization data-sets
% Submitted to Nature Methods in current form on 3.6.15, Manuscript number: NMETH-BC23142
% Authors: Matthew B. Stone* and Sarah L. Veatch+*
% * Department of Biophysics, University of Michigan, 930 N University, Ann Arbor MI 48109 (USA)
% + Corresponding Author: sveatch@umich.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dr =dr/calib; % convert bin size to pixel
max_r = 15*dr; %the maximum length to correlate out to 
rad_dist_bins = 0:dr:max_r; % the 1D spatial bins to histogram radii
Nframes = numel(data1);%getting the number of frames by the number of elements in the data1
%structure
mask_area = sum(sum(mask));
siz = size(mask); 
k = [1 cumprod(siz(1:end-1))]; %initializing for sub2ind workaround to 
%increase performance

all_dts = 1:max_tau+1; % real taus are 0:max_tau, to designate taus as 
% indexes in structures here using 1+taus to avoid indexing into 0

for t = all_dts 
    %initializing M(r) histograms at each time shift
    all_taus(t).M_11 = zeros(1,numel(rad_dist_bins));
    all_taus(t).M_22 = zeros(1,numel(rad_dist_bins));
    all_taus(t).M_12 = zeros(1,numel(rad_dist_bins));
    all_taus(t).M_21 = zeros(1,numel(rad_dist_bins));
    
    %zeroing N counter for each time shift
    all_taus(t).N_11 = 0;
    all_taus(t).N_22 = 0;
    all_taus(t).N_12 = 0;
    all_taus(t).N_21 = 0;
end
 

inner_rads = rad_dist_bins(1:end-1); % inner radius for each ring
outer_rads = rad_dist_bins(2:end); %outer radius for each ring
mean_rads = (inner_rads+outer_rads)/2; %center point of ring 
areas = pi*(outer_rads.^2 - inner_rads.^2); %area of each ring (in pixel^2)

%%

mask_rs = imresize(mask,1/dr); %resizing mask to calulate W(r)
% mask_rs has pixel size = dr, thus each bin of W_r = a correlation radii
I1 = double(mask_rs);         
rmax = max_r/dr;
L1 = size(I1, 1)+rmax*4; % size of fft2 (for zero padding)
L2 = size(I1, 2)+rmax*4; % size of fft2 (for zero padding)
w_r = real(fftshift(ifft2(abs(fft2(mask_rs, L1, L2)).^2))); % Normalization 
% for correct boundary conditions, w_r is the 2D autocorrelation of the mask

normalize = sum(sum(mask_rs)); %normalization factor for w_r
w_r_norm = w_r/normalize; %normalized 2D autocorrelation of mask
[W_r]  = radially_average(w_r_norm, rmax); %radially averaged 1D autocorrelation

%%

for j = 1:Nframes % going through all frames in user supplied data
    
    max_tau_j = min([max_tau Nframes-j]); % when within max_tau of last frame, 
    % this is the maximum tau that can be used
    
    dts = 0:max_tau_j; % actual tau values
    dt_inds = dts + 1; % tau values for indexes
    
    %applying mask to data from frame j
    [m_data1] = apply_mask([data1(j).x data1(j).y],mask,k);        
    [m_data2] = apply_mask([data2(j).x data2(j).y],mask,k);
    
    for t = dt_inds % going through all time lags (taus) with respect to frame j
    
        tshift = dts(t); % actual time lag to apply with respsect to frame j
        
        % data from frames thift after frame j
        DATA1_t = [data1(j+tshift).x data1(j+tshift).y];        
        DATA2_t = [data2(j+tshift).x data2(j+tshift).y];

        %applying mask to data from frame j+tshift
        [m_data1_t] = apply_mask(DATA1_t,mask,k);        
        [m_data2_t] = apply_mask(DATA2_t,mask,k);

        %computing the histogram M(r) and total number of radii N for channel 1
        %autocorrelation
        if ~isempty(m_data1) && ~isempty(m_data1_t) % if there are points in the mask                      
            num_combinations = size(m_data1_t,1)*size(m_data1,1); %total number of radii
            % between the two frames being correlated
            
           all_taus(t).N_11 = all_taus(t).N_11 + num_combinations;
           % keeping track of N, the total number of radii at each tau           
           
            [all_taus(t).M_11] = tabulate_radii(m_data1,m_data1_t,rad_dist_bins,all_taus(t).M_11);  
            % binning the radii between the two frames being correlated and
            % keeping track of the histogram of radii M(r) at each time shift                        
        end   
        
        %computing the histogram M and total number of radii N for channel 2
        %autocorrelation
        if ~isempty(m_data2) && ~isempty(m_data2_t)                           
            num_combinations = size(m_data2_t,1)*size(m_data2,1);
            all_taus(t).N_22 = all_taus(t).N_22 + num_combinations;
            [all_taus(t).M_22] = tabulate_radii(m_data2,m_data2_t,rad_dist_bins,all_taus(t).M_22);               
        end    
        
        %computing the histogram M and total number of radii N for channel
        %1->2 cross-correlation
        if ~isempty(m_data1) && ~isempty(m_data2_t)                          
            num_combinations = size(m_data1,1)*size(m_data2_t,1);
            all_taus(t).N_12 = all_taus(t).N_12 + num_combinations;
            [all_taus(t).M_12] = tabulate_radii(m_data1,m_data2_t,rad_dist_bins,all_taus(t).M_12);               
        end   
        
        %computing the histogram M and total number of radii N for channel
        %2->1 cross-correlation
        if ~isempty(m_data2) && ~isempty(m_data1_t)                          
            num_combinations = size(m_data2,1)*size(m_data1_t,1);  
            all_taus(t).N_21 = all_taus(t).N_21 + num_combinations;
            [all_taus(t).M_21] = tabulate_radii(m_data2,m_data1_t,rad_dist_bins,all_taus(t).M_21);               
        end                                
    end
end
%%

colorz = jet(numel(all_dts));

figure

for t = all_dts

    if t~=1 %not plotting autocorrelation at 0 time shift
    subplot(2,2,1)    
    C_11(t,:) = (all_taus(t).M_11(1:end-1)./(areas)./(W_r(1:numel(areas)))'*mask_area/(all_taus(t).N_11));
    dC_11(t,:) = sqrt((mask_area./areas./W_r(1:numel(areas))').^2.*(all_taus(t).M_11(1:end-1)./all_taus(t).N_11.^2).*(1+all_taus(t).M_11(1:end-1)./all_taus(t).N_11));     
    errorbar(mean_rads*calib,C_11(t,:),dC_11(t,:),'s-','MarkerFaceColor',colorz(t,:),'MarkerEdgeColor',colorz(t,:),'Color',colorz(t,:),'MarkerSize',2)
    hold on

    subplot(2,2,2)
    C_22(t,:) = (all_taus(t).M_22(1:end-1)./(areas)./(W_r(1:numel(areas)))'*mask_area/(all_taus(t).N_22));
    dC_22(t,:) = sqrt((mask_area./areas./W_r(1:numel(areas))').^2.*(all_taus(t).M_22(1:end-1)./all_taus(t).N_22.^2).*(1+all_taus(t).M_22(1:end-1)./all_taus(t).N_22));
    errorbar(mean_rads*calib,C_22(t,:),dC_22(t,:),'s-','MarkerFaceColor',colorz(t,:),'MarkerEdgeColor',colorz(t,:),'Color',colorz(t,:),'MarkerSize',2)
    hold on
    end
    
    subplot(2,2,3)
    C_12(t,:) = (all_taus(t).M_12(1:end-1)./(areas)./(W_r(1:numel(areas)))'*mask_area/(all_taus(t).N_12));
    dC_12(t,:) = sqrt((mask_area./areas./W_r(1:numel(areas))').^2.*(all_taus(t).M_12(1:end-1)./all_taus(t).N_12.^2).*(1+all_taus(t).M_12(1:end-1)./all_taus(t).N_12));
    errorbar(mean_rads*calib,C_12(t,:),dC_12(t,:),'s-','MarkerFaceColor',colorz(t,:),'MarkerEdgeColor',colorz(t,:),'Color',colorz(t,:),'MarkerSize',2)
    hold on
    
    subplot(2,2,4)
    C_21(t,:) = (all_taus(t).M_21(1:end-1)./(areas)./(W_r(1:numel(areas)))'*mask_area/(all_taus(t).N_21));
    dC_21(t,:) = sqrt((mask_area./areas./W_r(1:numel(areas))').^2.*(all_taus(t).M_21(1:end-1)./all_taus(t).N_21.^2).*(1+all_taus(t).M_21(1:end-1)./all_taus(t).N_21));
    errorbar(mean_rads*calib,C_21(t,:),dC_21(t,:),'s-','MarkerFaceColor',colorz(t,:),'MarkerEdgeColor',colorz(t,:),'Color',colorz(t,:),'MarkerSize',2)
    hold on


end

subplot(2,2,1)  
ylabel('G(r,\tau)')
xlabel('interparticle distance (nm)')
title('data1 autocorrelation','fontweight','bold')

subplot(2,2,2)  
ylabel('G(r,\tau)')
xlabel('interparticle distance (nm)')
title('data2 autocorrelation','fontweight','bold')

subplot(2,2,3)  
ylabel('C(r,\tau)')
xlabel('interparticle distance (nm)')
title('data1(t) x data2(t+\tau)','fontweight','bold')

subplot(2,2,4)  
ylabel('C(r,\tau)')
xlabel('interparticle distance (nm)')
title('data2(t) x data1(t+\tau)','fontweight','bold')

if numel(all_dts)>1 % if there is autocorrelaiton data
    corrs.C_11 = C_11; 
    corrs.dC_11 = dC_11; 

    corrs.C_22 = C_22; 
    corrs.dC_22 = dC_22; 
end
    
corrs.C_12 = C_12; 
corrs.dC_12 = dC_12; 

corrs.C_21 = C_21; 
corrs.dC_21 = dC_21; 

params.radii = mean_rads*calib; 
end %returns OUTPUT


function [m_data1] = apply_mask(data1,mask,k)

if isempty(data1); %returns empty if no data 
    m_data1 = [];
return
end

X1 = data1(:,1); % x positions of data
Y1 = data1(:,2); % y positions of data

inders = 1:numel(X1); %numbered indexes of original data

NN = X1<size(mask,1); % getting rid of out of range indexes
WW = X1>0;
NM = Y1<size(mask,2); 
WA = Y1>0;

RERE = NN.*WW.*NM.*WA;
loggg1 = logical(RERE); %logical indexes of data range of mask
        
if sum(loggg1) == 0; %if there is no data, return everything empty
    m_data1 = [];     
return
end        
            
subs = [ceil(double(X1(loggg1))) ceil(double(Y1(loggg1)))];
% the x and y positions that are within the range of the mask
% subs has removed particles that are out of range of the mask      
% subs is ceiling in order to index directly into the mask 

struct_inds = inders(loggg1);
%the numerical indexes of the original data in range of mask

ch1_inds = 1 + (subs(:,1)-1)*k(1) + (subs(:,2)-1)*k(2);
% getting 1D index from 2D subscripts subs [x y] using a 
% workaround for sub2ind. sub2ind calls led to decreased performance        

keep1 = mask(ch1_inds); % applying the logical mask to the data     
inds_zeros = struct_inds'.*keep1; % indexes of all data within mask
% with zeros at indexes of data outside of mask        

final_inds = inds_zeros(inds_zeros>0);        
m_data1 = data1(final_inds,:); %the masked data1                      
end

function [vals] = radially_average(I, rmax)
if nargin<2, rmax = 100; end

% finds the center of the image I in each dimension
center1 = ceil((size(I, 1)+1)/2); 
center2 = ceil((size(I, 2)+1)/2);

% sets the interval to calculate the radial average, 
% centered on image I center pixel
range1 = center1-rmax:center1+rmax;
range2 = center2-rmax:center2+rmax;

% the values from the interval
zvals = I(range1, range2);

% creating a meshgrid from the interval
[xvals yvals] = meshgrid([-rmax:rmax],[-rmax:rmax]);

% transform to polar coordinates with v as image values 
[theta,r,v] = cart2pol(xvals,yvals, zvals);

% arrange from small to large r values
Ar = reshape(r,1, numel(r));
[rr,ind] = sort(Ar);

% the corresponding image values
Avals = reshape(v,1, numel(v));
vv = Avals(ind);

% bin by radii and average values in bins
r = 0:floor(max(rr));
[n, bin] = histc(rr, r-.5);
vals = zeros(rmax+1,1); 
for j = 1:rmax+1;%length(r),
    m = bin==j;
    n2 = sum(m);
    if n2==0, vals(j)=0; er(j)=0; 
    else
        vals(j) = sum(m.*vv)/n2;        
    end
end

end

function [histobam] = tabulate_radii(DATA1,DATA2,rad_dist_bins,histobam)

m_data1 = DATA1; %the masked data1
m_data2 = DATA2; %the masked data2  



for k = 1: numel(m_data1(:,1)) %cycling through only data points within mask

        A = m_data1(k,:); %particle k's x and y position

        B = A(ones(size(m_data2,1),1),:); %repeat particle 1's position 
        %to match size of data2
            diff = B - m_data2; %comparing particle k to all particles 
            %in data2
            diff_squared = diff.*diff; % [x_diff^2 y_diff^2]
            r_diffs = (sqrt(diff_squared(:,1) + diff_squared(:,2)))'; 
            %radial difference values



[new_hists indexers2] = histc(r_diffs,rad_dist_bins); %binning by radius 
% rad_dist_bins = radial distribution 

%indexing the column indexes by k, (m_data1(k))
%amazingly, this index takes car of both channels indexing in an elegant
%way: all_radii_indexes(k,j) will return a vector containing which 
%bins the radii between m_data(k) and m_data(j) ((NOTE m_data is actually 
% N by 2 since it contains x and y values for position...)

histobam = histobam + new_hists; %adding the histograms for g(r)
  
end

end
