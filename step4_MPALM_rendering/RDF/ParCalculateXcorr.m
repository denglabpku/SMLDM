function [corrs, params] = ParCalculateXcorr(data1_x, data1_y, data2_x, data2_y, mask, calib, dr, num_pts, item)
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
% NUM_PTS is the MULTIPLY for maximum length to correlate out to 
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
% Function requirement: radially_average.m, apply_mask.m, tabulate_radii.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 9
    item = [1, 1, 1];
end

data1(1).x = data1_x;
data1(1).y = data1_y;
data2(1).x = data2_x;
data2(1).y = data2_y;

dr =dr/calib; % convert bin size to pixel
max_r = num_pts * dr; % max_r is the maximum length to correlate out to 
rad_dist_bins = 0:dr:max_r; % the 1D spatial bins to histogram radii
Nframes = numel(data1); % getting the number of frames by the number of elements in the data1, always when pool locs from all frames
%structure
mask_area = sum(sum(mask));
siz = size(mask);
k = [1 cumprod(siz(1:end-1))]; %initializing for sub2ind workaround to 
%increase performance

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
    %initializing M(r) histograms at each time shift
    all_corrs(j).M_11 = zeros(1,numel(rad_dist_bins));
    all_corrs(j).M_22 = zeros(1,numel(rad_dist_bins));
    all_corrs(j).M_12 = zeros(1,numel(rad_dist_bins));
    
    %zeroing N counter for each time shift
    all_corrs(j).N_11 = 0;
    all_corrs(j).N_22 = 0;
    all_corrs(j).N_12 = 0;

    %applying mask to data from frame j
    [m_data1] = apply_mask([data1(j).x data1(j).y],mask,k);        
    [m_data2] = apply_mask([data2(j).x data2(j).y],mask,k);

    %computing the histogram M(r) and total number of radii N for channel 1
    % and channel 2
    % autocorrelation and cross-correlation 1<->2
    if ~isempty(m_data1) && ~isempty(m_data2)
        if item(1) > 0
            num_combinations = size(m_data1,1)^2;
            all_corrs(j).N_11 = all_corrs(j).N_11 + num_combinations;
            [all_corrs(j).M_11] = tabulate_radii(m_data1,m_data1,rad_dist_bins,all_corrs(j).M_11); 
        end
        if item(2) > 0
            num_combinations = size(m_data2, 1)^2;
            all_corrs(j).N_22 = all_corrs(j).N_22 + num_combinations;
            [all_corrs(j).M_22] = tabulate_radii(m_data2,m_data2,rad_dist_bins,all_corrs(j).M_22);
        end
        if item(3) > 0
            num_combinations = size(m_data1, 1) * size(m_data2, 1);
            all_corrs(j).N_12 = all_corrs(j).N_12 + num_combinations;
            [all_corrs(j).M_12] = tabulate_radii(m_data1,m_data2,rad_dist_bins,all_corrs(j).M_12);
        end
    end 

    if (item(1) > 0) && (j < Nframes)
        all_corrs(j).t_11 = zeros(1,numel(rad_dist_bins));
        all_corrs(j).n_11 = 0;
        [m_data1_t] = apply_mask([data1(j+1).x data1(j+1).y],mask,k);

        num_combinations = size(m_data1, 1) * size(m_data1_t, 1);
        all_corrs(j).n_11 = all_corrs(j).n_11 + num_combinations;
        [all_corrs(j).t_11] = tabulate_radii(m_data1,m_data1_t,rad_dist_bins,all_corrs(j).t_11);
    end
end
%%
colorz = jet(Nframes);

figure

for t = 1:Nframes

    if (item(1) > 0)
        subplot(2,2,1)    
        C_11(t,:) = (all_corrs(t).M_11(1:end-1)./(areas)./(W_r(1:numel(areas)))'*mask_area/(all_corrs(t).N_11));
        dC_11(t,:) = sqrt((mask_area./areas./W_r(1:numel(areas))').^2.*(all_corrs(t).M_11(1:end-1)./all_corrs(t).N_11.^2).*(1+all_corrs(t).M_11(1:end-1)./all_corrs(t).N_11));     
        errorbar(mean_rads*calib,C_11(t,:),dC_11(t,:),'s-','MarkerFaceColor',colorz(t,:),'MarkerEdgeColor',colorz(t,:),'Color',colorz(t,:),'MarkerSize',2)
        hold on
        if t < Nframes
            subplot(2,2,4)
            T_11(t,:) = (all_corrs(t).t_11(1:end-1)./(areas)./(W_r(1:numel(areas)))'*mask_area/(all_corrs(t).n_11));
            dT_11(t,:) = sqrt((mask_area./areas./W_r(1:numel(areas))').^2.*(all_corrs(t).t_11(1:end-1)./all_corrs(t).n_11.^2).*(1+all_corrs(t).t_11(1:end-1)./all_corrs(t).n_11));
            errorbar(mean_rads*calib,T_11(t,:),dT_11(t,:),'s-','MarkerFaceColor',colorz(t,:),'MarkerEdgeColor',colorz(t,:),'Color',colorz(t,:),'MarkerSize',2)
            hold on
        end
    end

    if (item(2) > 0)
        subplot(2,2,2)
        C_22(t,:) = (all_corrs(t).M_22(1:end-1)./(areas)./(W_r(1:numel(areas)))'*mask_area/(all_corrs(t).N_22));
        dC_22(t,:) = sqrt((mask_area./areas./W_r(1:numel(areas))').^2.*(all_corrs(t).M_22(1:end-1)./all_corrs(t).N_22.^2).*(1+all_corrs(t).M_22(1:end-1)./all_corrs(t).N_22));
        errorbar(mean_rads*calib,C_22(t,:),dC_22(t,:),'s-','MarkerFaceColor',colorz(t,:),'MarkerEdgeColor',colorz(t,:),'Color',colorz(t,:),'MarkerSize',2)
        hold on
    end
    
    if (item(3) > 0)
        subplot(2,2,3)
        C_12(t,:) = (all_corrs(t).M_12(1:end-1)./(areas)./(W_r(1:numel(areas)))'*mask_area/(all_corrs(t).N_12));
        dC_12(t,:) = sqrt((mask_area./areas./W_r(1:numel(areas))').^2.*(all_corrs(t).M_12(1:end-1)./all_corrs(t).N_12.^2).*(1+all_corrs(t).M_12(1:end-1)./all_corrs(t).N_12));
        errorbar(mean_rads*calib,C_12(t,:),dC_12(t,:),'s-','MarkerFaceColor',colorz(t,:),'MarkerEdgeColor',colorz(t,:),'Color',colorz(t,:),'MarkerSize',2)
        hold on
    end

end

if item(1) > 0
    subplot(2,2,1)  
    ylabel('G(r,\tau)')
    xlabel('interparticle distance (nm)')
    title('data1 autocorrelation','fontweight','bold')
    corrs.C_11 = C_11; 
    corrs.dC_11 = dC_11; 
    
    subplot(2,2,4)  
    ylabel('C(r,\tau)')
    xlabel('interparticle distance (nm)')
    title('data1(t) x data1(t+1)','fontweight','bold')
    
    if Nframes > 1
        corrs.T_11 = T_11; 
        corrs.dT_11 = dT_11; 
    end
end

if item(2) > 0
    subplot(2,2,2)  
    ylabel('G(r,\tau)')
    xlabel('interparticle distance (nm)')
    title('data2 autocorrelation','fontweight','bold')
    corrs.C_22 = C_22; 
    corrs.dC_22 = dC_22; 
end

if item(3) > 0
    subplot(2,2,3)  
    ylabel('C(r,\tau)')
    xlabel('interparticle distance (nm)')
    title('data1 x data2','fontweight','bold')
    corrs.C_12 = C_12; 
    corrs.dC_12 = dC_12; 
end

params.radii = mean_rads*calib; 
end %returns OUTPUT