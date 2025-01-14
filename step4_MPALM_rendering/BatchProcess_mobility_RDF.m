% This script is to analyze MPALM mobility distribution plus membrane
% protein RDF.

% The input file is from BatchInput_mobility_RDF.m

%% Fitting averaged beads image with Gaussian function
% sigma is used in motion-blur simulation
close all
root_dir = 'E:\PROCESS-SPT\2023\simPSF_results\DeepSTORMDataset\20240321_debug2_shotnoise';
img_file = dir(fullfile(root_dir,'*.tif'));
name_contains ="TIRF_average_psf_qd655"; 
pixel_size = 0.11; % unit um
valid_file = contains({img_file.name},name_contains);
val_img_file = img_file(valid_file);

fitresult = nan(1,5);


for i = 1:1
    img = imread(fullfile(val_img_file.folder,val_img_file.name),i);
    [imHeight,imWidth] = size(img);

    % Assuming your image is stored in the variable 'image_data'
    % Convert the 2D image into separate X and Y coordinate matrices
    [X, Y] = meshgrid(1:imHeight, 1:imWidth);
    
    % Reshape the X, Y, and image data to be compatible with prepareSurfaceData
    X = reshape(X, [], 1);
    Y = reshape(Y, [], 1);
    intensity_values = double(reshape(img, [], 1));
    
    % Now use prepareSurfaceData
    [Xq, Yq, Zq] = prepareSurfaceData(X, Y, intensity_values);
    
    fitresult(i,:) = cell_fit_gaussian2(Xq,Yq,Zq,pixel_size);

    xx = linspace(1,imHeight,100);
    f = @(I,m,r0,x0,x) (I * exp(-(1/(2*r0.^2))*(x - x0).^2) + m);
    yy = f(fitresult(i,1),fitresult(i,2),fitresult(i,3),fitresult(i,4),xx);

    fig = figure("Position",[202 351 918 270]); tiledlayout('flow','Padding','tight')
    nexttile;
    imshow(img,[],'InitialMagnification',1000)
    title("Averaged PSF at focal plane")

    nexttile;hold on;
    plot((1:imHeight)*pixel_size,img(:,17),'--s'); 
    plot(xx*pixel_size,yy,'r-')
    legend(["signal","Gauss fit"])
    title(sprintf("Profile along y axis, r0 : %.4f",fitresult(i,3)))

    nexttile;hold on;
    plot((1:imHeight)*pixel_size,img(17,:),'--s'); 
    plot(xx*pixel_size,yy,'r-')
    legend(["signal","Gauss fit"])
    title(sprintf("Profile along x axis, r0 : %.4f",fitresult(i,3)))

end
% exportgraphics(fig,"beads_fitting.pdf")

%% Prepare the multi-normal fitting model using simulated dataset
% ---------- old simulation results --------------%
% tab_psf_fitresult = readtable('E:\PROCESS-SPT\2023\simPSF_results\UNetDetectionDataset\20231204_NoNoiseNoLocError_30ms_160nmPixSize_5k\NoNoise_SmoothArea_PixelArea.csv');
% tab_psf_fitresult = readtable('E:\PROCESS-SPT\2023\simPSF_results\UNetDetectionDataset\20231207_PSNR25_5kImgs_110nmPixSize_beadsPSF\tab_psf_fitresult.csv');
% tab_psf_fitresult = readtable('E:\PROCESS-SPT\2023\simPSF_results\UNetDetectionDataset\20231207_PSNR25_5kImgs_110nmPixSize\tab_psf_fitresult.csv');
% tab_psf_fitresult = readtable('E:\PROCESS-SPT\2023\simPSF_results\UNetDetectionDataset\20231207_NoNoiseNoLocError_30ms_110nmPixSize_5k_beadsPSF\NoNoise_SmoothArea_PixelArea.csv'); 
% tab_psf_fitresult = readtable('E:\PROCESS-SPT\2023\simPSF_results\UNetDetectionDataset\20240111_SNR19-35_NoLocError_110nmPixSize\MSD.csv');
% tab_psf_fitresult = readtable('E:\PROCESS-SPT\2023\simPSF_results\UNetDetectionDataset\20240111_SNR19-35_NoLocError_110nmPixSize_log10D\MSD.csv');
% tab_psf_fitresult = readtable('D:\PROCESS-SPT\2023\simPSF_results\UNetDetectionDataset\20240219_SNR19-35_NoLocError_110nmPixSize_log10D\MSD.csv'); % more data (1000 per D) in small D
% tab_psf_fitresult = readtable('E:\PROCESS-SPT\2023\simPSF_results\DeepSTORMDataset\20240321_debug2_shotnoise\MBX_20231220_110nmPix_rep2_epoch9_r1p42_SRArea.csv'); % debug

% D = [0.00100    0.00158	0.00251	0.00398	0.00631	0.01000	0.01585	0.02512	0.03981	0.06310 0.10000];
% D = [0.00100    0.00158	0.00251	0.00398	0.00631	0.01000	0.01585	0.02512	0.03981	0.06310	0.10000	0.15849	0.25119	0.39811	0.63096	1.00000	1.58489	2.51189	3.98107	6.30957	7.94328	10.00000	12.58925	15.84893	19.95262	25.11886	31.62278];
% D = [0.01000	0.01585	0.02512	0.03981	0.06310	0.10000	0.15849	0.25119	0.39811	0.63096	1.00000	1.58489	2.51189	3.98107	6.30957	7.94328	10.00000	12.58925	15.84893	19.95262	25.11886	31.62278];
% D = [0.10000	0.15849	0.25119	0.39811	0.63096	1.00000	1.58489	2.51189	3.98107	6.30957	7.94328	10.00000	12.58925	15.84893	19.95262	25.11886	31.62278];

clear;close all; clc;
% ---------- latest simulation results --------%
tab_psf_fitresult = readtable('E:\PROCESS-SPT\2023\simPSF_results\DeepSTORMDataset\20240321_debug2_shotnoise\MBX_20231220_110nmPix_rep2_epoch9_r1p41long_SRArea.csv'); % latest simulation

D = [0.0666 0.10000	0.15849	0.25119	0.39811	0.63096	1.00000	1.58489	2.51189	3.98107	6.30957	7.94328	10.00000	12.58925	15.84893	19.95262	25.11886	31.62278];
mu = zeros(size(D));
sigma = zeros(size(D));
mu(1) = 2.1; % for beads fitting

% % discard D below 0.1
% tab_psf_fitresult = tab_psf_fitresult(tab_psf_fitresult.DiffCoeff >= 0.10000,:);

% pool all results below 0.1, which cannot distinguish
tab_psf_fitresult.DiffCoeff(tab_psf_fitresult.DiffCoeff < 0.10000,:) = 0.0666;

% Visulize multiple Gaussian distribution
fig1 = figure('Position',[22 580 800 314]); tiledlayout(1,2);
color_pallete = colormap("jet");
color_idx = round(1:length(color_pallete)/length(D):length(color_pallete));
nexttile;hold on;
for i = 1:length(D)      
    % Estimate the PDF using kernel smoothing
    [f, xi] = ksdensity(log10(tab_psf_fitresult.mask_areas1(tab_psf_fitresult.DiffCoeff == D(i))),'Function','CDF');
    % Plot the kernel smoothed probability
    plot(xi,f,'Color',color_pallete(color_idx(i),:),'LineWidth',0.5)
end
hold off
xlabel('log10(SR Area) (unit: pix)')
ylabel('CDF')
title('Simulated area')
xlim([1.5 4.2])

nexttile; hold on;
bandwidth = [0.1*ones(1,14),0.05*ones(1,4)];
for i = 1:length(D)      
    % Estimate the PDF using kernel smoothing
    [f, xi] = ksdensity(log10(tab_psf_fitresult.mask_areas1(tab_psf_fitresult.DiffCoeff == D(i))),'Bandwidth',bandwidth(i));
    % Plot the kernel smoothed probability
    plot(xi,f,'Color',color_pallete(color_idx(i),:),'LineWidth',2)
    pd = fit(xi',f','gauss1');
    mu(i) = pd.b1;
    sigma(i) = pd.c1;
   
end
ldg = legend(string(D));
lgd.Layout.Tile = 'east';
hold off
xlabel('log10(SR Area) (unit: pix)')
ylabel('PDF')
title('Simulated area')
xlim([1.5 4.2])
% exportgraphics(fig,fullfile(table_path,'PDF_AveDist_manyD_at30ms_SNR22.pdf'))

% Visualize the raw sim data
fig2 = figure('Position',[205 432 1007 304]); tiledlayout(3,6);
color_pallete = colormap("jet");
color_idx = round(1:length(color_pallete)/length(D):length(color_pallete));
for i = 1:length(D)      
    nexttile;
    f = log10(tab_psf_fitresult.mask_areas1(tab_psf_fitresult.DiffCoeff == D(i)));
    histogram(f,linspace(1,4,20),'EdgeColor',color_pallete(color_idx(i),:),'FaceColor','none','Normalization','probability')
    title(sprintf('D: %.3f',D(i)));
    xlim([1.5 4.5])
end

%% Set up the fitting model (for multi-state Gaussian model only)
fprintf('Loading fitting model ...\n')
% sigma is a fitting param
cdf_pd = @(i,sigma,x) cdf('Normal',x,mu(i),sigma);
ft_cdf = @(w,x) ...
    w(1) * cdf_pd(1,w(19),x) + w(2) * cdf_pd(2,w(20),x) + ...
    w(3) * cdf_pd(3,w(21),x) + w(4) * cdf_pd(4,w(22),x) + ...
    w(5) * cdf_pd(5,w(23),x) + w(6) * cdf_pd(6,w(24),x) + ...
    w(7) * cdf_pd(7,w(25),x) + w(8) * cdf_pd(8,w(26),x) + ...
    w(9) * cdf_pd(9,w(27),x) + w(10) * cdf_pd(10,w(28),x) + ...
    w(11) * cdf_pd(11,w(29),x) + w(12) * cdf_pd(12,w(30),x) + ...
    w(13) * cdf_pd(13,w(31),x) + w(14) * cdf_pd(14,w(32),x) + ...
    w(15) * cdf_pd(15,w(33),x) + w(16) * cdf_pd(16,w(34),x) + ...
    w(17) * cdf_pd(17,w(35),x) + w(18) * cdf_pd(18,w(36),x);

pdf_pd = @(i,sigma,x) pdf('Normal',x,mu(i),sigma);
ft_pdf = @(w,x) ...
    w(1) * pdf_pd(1,w(19),x) + w(2) * pdf_pd(2,w(20),x) + ...
    w(3) * pdf_pd(3,w(21),x) + w(4) * pdf_pd(4,w(22),x) + ...
    w(5) * pdf_pd(5,w(23),x) + w(6) * pdf_pd(6,w(24),x) + ...
    w(7) * pdf_pd(7,w(25),x) + w(8) * pdf_pd(8,w(26),x) + ...
    w(9) * pdf_pd(9,w(27),x) + w(10) * pdf_pd(10,w(28),x) + ...
    w(11) * pdf_pd(11,w(29),x) + w(12) * pdf_pd(12,w(30),x) + ...
    w(13) * pdf_pd(13,w(31),x) + w(14) * pdf_pd(14,w(32),x) + ...
    w(15) * pdf_pd(15,w(33),x) + w(16) * pdf_pd(16,w(34),x) + ...
    w(17) * pdf_pd(17,w(35),x) + w(18) * pdf_pd(18,w(36),x);

% % Set up fittype and options.
% pdf_pd = @(i,x) pdf('Normal',x,mu(i),sigma(i));
% ft_pdf = @(w,x) ...
%     w(1) * pdf_pd(1,x) + w(2) * pdf_pd(2,x) + ...
%     w(3) * pdf_pd(3,x) + w(4) * pdf_pd(4,x) + ...
%     w(5) * pdf_pd(5,x) + w(6) * pdf_pd(6,x) + ...
%     w(7) * pdf_pd(7,x) + w(8) * pdf_pd(8,x) + ...
%     w(9) * pdf_pd(9,x) + w(10) * pdf_pd(10,x) + ...
%     w(11) * pdf_pd(11,x) + w(12) * pdf_pd(12,x) + ...
%     w(13) * pdf_pd(13,x) + w(14) * pdf_pd(14,x) + ...
%     w(15) * pdf_pd(15,x) + w(16) * pdf_pd(16,x) + ...
%     w(17) * pdf_pd(17,x) + w(18) * pdf_pd(18,x);
% 
% 
% % Set up fittype and options.
% cdf_pd = @(i,x) cdf('Normal',x,mu(i),sigma(i));
% ft_cdf = @(w,x) ...
%     w(1) * cdf_pd(1,x) + w(2) * cdf_pd(2,x) + ...
%     w(3) * cdf_pd(3,x) + w(4) * cdf_pd(4,x) + ...
%     w(5) * cdf_pd(5,x) + w(6) * cdf_pd(6,x) + ...
%     w(7) * cdf_pd(7,x) + w(8) * cdf_pd(8,x) + ...
%     w(9) * cdf_pd(9,x) + w(10) * cdf_pd(10,x) + ...
%     w(11) * cdf_pd(11,x) + w(12) * cdf_pd(12,x) + ...
%     w(13) * cdf_pd(13,x) + w(14) * cdf_pd(14,x) + ...
%     w(15) * cdf_pd(15,x) + w(16) * cdf_pd(16,x) + ...
%     w(17) * cdf_pd(17,x) + w(18) * cdf_pd(18,x);

%% Run mobility distribution fitting
close all;
D_SR_const = 1162 ;%2007; % 1162

draw_roi = false;
load_roi = false;

% loop through each sample
for iSamp = 1:length(data_struct)
    fprintf('Loading sample: %s \n',data_struct(iSamp).SampleName);
    tab_all_SR_pred = table();
    temp_single_data = {};
    tot_cell_num = length(data_struct(iSamp).workspaces);

    % loop through each cell
    for iCell = 1:tot_cell_num
        clear tab_SR_pred
        % SRFile = dir(fullfile(data_struct(iSamp).workspaces(iCell).folder,data_struct(iSamp).workspaces(iCell).name, 'UNet_*.csv'));
        % tab_SR_pred = readtable(fullfile(SRFile.folder,SRFile.name));
        tab_SR_pred = readtable(fullfile(data_struct(iSamp).workspaces(iCell).folder,data_struct(iSamp).workspaces(iCell).name, data_struct(iSamp).SRFile));
        if draw_roi
            fig = figure('Position',[471 104 822 738]);
            scatter(tab_SR_pred.Xpos, tab_SR_pred.Ypos,2,'filled','MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none','MarkerFaceAlpha',0.3);
            xlabel('X Position');
            ylabel('Y Position');
            title('Select ROI');
            axis image ij
            
            % Draw a polygon ROI
            hROI = drawpolygon('LineWidth', 2, 'Color', 'r');
            
            % Get the X and Y coordinates of the vertices of the polygon
            roiVertices = hROI.Position;
            
            % Use inpolygon to find points within the ROI
            [in, on] = inpolygon(tab_SR_pred.Xpos, tab_SR_pred.Ypos, roiVertices(:,1), roiVertices(:,2));
            
            % Filter the table to include only rows within the ROI
            tab_SR_pred = tab_SR_pred(in | on, :);

            % save(fullfile(fullfile(data_struct(iSamp).workspaces(iCell).folder,data_struct(iSamp).workspaces(iCell).name, ...
            %     sprintf('roi_metadata.mat'))),'hROI');

            save(fullfile(SRFile.folder,sprintf('roi_metadata.mat')),'hROI');
            close(fig)
        end

        if load_roi

            load(fullfile(SRFile.folder,sprintf('roi_metadata.mat')),'hROI');

            % Get the X and Y coordinates of the vertices of the polygon
            roiVertices = hROI.Position;

            fig = figure('Position', [471 104 822 738],'Visible','off');            
            scatter(tab_SR_pred.Xpos, tab_SR_pred.Ypos, 2, 'filled', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.3);
            xlabel('X Position');
            ylabel('Y Position');
            title('Select ROI');
            axis image ij
    
            % Draw a polygon ROI
            hROI = drawpolygon('Position',roiVertices,'LineWidth', 2, 'Color', 'r');
                    
            % Use inpolygon to find points within the ROI
            [in, on] = inpolygon(tab_SR_pred.Xpos, tab_SR_pred.Ypos, roiVertices(:, 1), roiVertices(:, 2));
    
            % Filter the table to include only rows within the ROI
            tab_SR_pred = tab_SR_pred(in | on, :);

            
            close(fig)
        end

        temp_single_data{iCell} = tab_SR_pred.SRArea;
        if ~isempty(find(tab_SR_pred.SRArea == -Inf, 1))
            uiwait;
        end
        
        fprintf('localization number %d ...\n',height(tab_SR_pred))
        tab_all_SR_pred = [tab_all_SR_pred;tab_SR_pred];
        FinalResults{iSamp}.merge_data = tab_all_SR_pred.SRArea;
        FinalResults{iSamp}.merge_D = (tab_all_SR_pred.SRArea./D_SR_const).^2;
        % FinalResults{iSamp}.locPerFrame = [FinalResults{iSamp}.locPerFrame height(tab_SR_pred)/max(tab_SR_pred.Frame)];

        if FitParams.DoSingleFit
            % single cell multi Burr fitting            
            single_D = (temp_single_data{iCell}./D_SR_const).^2;
            single_D = single_D(single_D>0); % remove error estimation
            fprintf("Fitting dataset: %s, Localization number %d ...\n",data_struct(iSamp).workspaces(iCell).name,length(single_D));
            [FinalResults{iSamp}.single_fitResult{iCell},~,FinalResults{iSamp}.single_fitOutput{iCell}] = two_gauss_fit(log10(single_D));

            % normalize to sum 1 for each fraction
            w = coeffvalues(FinalResults{iSamp}.single_fitResult{iCell});
            p1 = w(1)/(w(1)+w(4));p2 = w(4)/(w(1)+w(4));  w(1) =p1;  w(4) = p2;            
            FinalResults{iSamp}.single_coefficient(iCell,:) = [p1 10^w(2) p2 10^w(5)]; % [component1 faction, comp1 mean, component2 faction, comp2 mean]

            % previous version
            % [FinalResults{iSamp}.single_w(iCell,:),~,FinalResults{iSamp}.residuals{iCell},~,~] = multi_gauss_fit(log10(temp_single_data{iCell}),ft_cdf);
            % FinalResults{iSamp}.single_norm_w(iCell,:) = FinalResults{iSamp}.single_w(iCell,1:18)/sum(FinalResults{iSamp}.single_w(iCell,1:18));
        end
     end
    
    if FitParams.DoMergeFit
        fprintf("Fitting merged dataset of above ...\n");
        % single cell multi Burr fitting            
        pooled_D = (FinalResults{iSamp}.merge_data./D_SR_const).^2;
        pooled_D = pooled_D(pooled_D>0); % remove error estimation
        [FinalResults{iSamp}.merge_fitResult,~,FinalResults{iSamp}.merge_fitOutput] = two_gauss_fit(log10(pooled_D));

        % normalize to sum 1 for each fraction
        w = coeffvalues(FinalResults{iSamp}.merge_fitResult);
        p1 = w(1)/(w(1)+w(4));p2 = w(4)/(w(1)+w(4));  w(1) =p1;  w(4) = p2;            
        FinalResults{iSamp}.merge_coefficient = [p1 10^w(2) p2 10^w(5)]; % [component1 faction, comp1 mean, component2 faction, comp2 mean]

        % previous version
        % all cell multi normal fitting
        % [FinalResults{iSamp}.w,~,~,~,~] = multi_gauss_fit(log10(FinalResults{iSamp}.merge_data),ft_cdf);
        % FinalResults{iSamp}.merge_norm_w = FinalResults{iSamp}.w(1:18);
        % FinalResults{iSamp}.merge_norm_w = FinalResults{iSamp}.merge_norm_w/sum(FinalResults{iSamp}.merge_norm_w);
    end
    
    % ------------ singe-cell results  ------------ % 
    if FitParams.DoSinglePlot
  
        % ===================== single cell PLOT ============================= %

        % PDF of indivisual cell
        figure('Position',[375 641 560 159]);
        color_pallete = colormap("jet");
        color_idx = round(1:length(color_pallete)/tot_cell_num:length(color_pallete));
        hold on;
        for iCell = 1:tot_cell_num 
%             xgrid = 0:200:10000;
%             pd = fitdist(temp_single_data{iCell},'Kernel');
%             y = pdf(pd,xgrid);
%             plot(xgrid,y,'LineWidth',2,'Color',color_pallete(color_idx(iCell),:))
            
            xgrid = linspace(0,4,100);
            h = histogram(log10(temp_single_data{iCell}),xgrid,'Normalization','PDF','Visible','off','HandleVisibility','off');
            xgrid = (xgrid(1:end-1)+xgrid(2:end))/2; 
            plot(xgrid,h.Values,'LineWidth',2,'Color',color_pallete(color_idx(iCell),:))
        end
        xlabel('log10(PT area) (unit: A.U.)')
        ylabel('PDF');xlim([1.5 4]);
        title('single cell overlay')
        
        if FitParams.DoSingleFit
            % plot individual cell fitting results
            counter = 0;  % Counter variable to keep track of the number of plots
            figure('Units','inches','Position',[2.6354 5.6667 15.0000 2.9271]);  % Create the first figure
            tiledlayout(2,2);
            for iCell = 1:tot_cell_num
                % Check if counter is a multiple of 8figure;
                if mod(counter, 4) == 0 && counter > 0
                    figure('Units','inches','Position',[2.6354 5.6667 15.0000 2.9271]);  % Create a new figure
                    tiledlayout(2,2);
                    counter = 0;
                end
                nexttile(counter+1);hold on;
                xgrid = -3:0.1:2;
                histogram(log10(single_D),xgrid,'Normalization','PDF')
                h = plot(FinalResults{iSamp}.single_fitResult{iCell});
                % xgrid = linspace(0,4,100);
                % histogram(log10(temp_single_data{iCell}),xgrid,'Normalization','PDF','FaceColor','none')
                % plot(xgrid,ft_pdf(FinalResults{iSamp}.single_w(iCell,:),xgrid),'LineWidth',2);
                xlabel( 'Log10(PT Area) (unit: a.u.)', 'Interpreter', 'none' );
                ylabel( 'PDF', 'Interpreter', 'none' );
                title(sprintf('%s',data_struct(iSamp).workspaces(iCell).name),'Interpreter','none');
                xlim([1.5 4])
                counter = counter + 1;  % Increment the counter
            end
    
            % plot individual cell fitting residuals
            counter = 0;  % Counter variable to keep track of the number of plots
            figure('Units','inches','Position',[2.6354 4.9062 15 3.6875]);  % Create the first figure
            tiledlayout(2,2);
            % plot individual cell fitting results
            for iCell = 1:tot_cell_num
                % Check if counter is a multiple of 8
                if mod(counter, 4) == 0 && counter > 0
                    figure('Units','inches','Position',[2.6354 4.9062 15 3.6875]);  % Create a new figure
                    tiledlayout(2,2);
                    counter = 0;
                end
    
                nexttile(counter+1);hold on;
                xgrid = linspace(0,4,100);
                h = plot(xgrid,FinalResults{iSamp}.residuals{iCell});
                h.LineWidth = 2;
                xlabel( 'PT Area', 'Interpreter', 'none' );
                ylabel( 'Residuals', 'Interpreter', 'none' );
                title(sprintf('%s',data_struct(iSamp).workspaces(iCell).name),'Interpreter','none');
                xlim([1.5 4])
    
                counter = counter + 1;  % Increment the counter
            end
    
            % heatmap of weights at diffusion coeffcients
            figure;
            h = heatmap(FinalResults{iSamp}.single_norm_w);
            h.XData = D; h.ColorScaling = 'scaledrows';
            colormap parula
            xlabel('Diffusion coefficient')
            title('Heatmap of normalized weights in each sample')
        end

%         % plot overlaid single cell weight
%         % Define colors for each curve
%         color_pallete = colormap("jet");
%         color_idx = round(1:length(color_pallete)/tot_cell_num:length(color_pallete));
%         figure; hold on
%         for iCell = 1:tot_cell_num
%             plot(D,FinalResults{iSamp}.single_norm_w(iCell,:),'LineWidth',2,'Color',color_pallete(color_idx(iCell),:))
%         end
%         ax = gca;
%         ax.XScale = 'log';
%         legend;
%         xlabel('Diffusion coefficient')
%         ylabel('Normalized weight')
%         title('single cell multi-normal fitting from motion blur')
        
%         % ===================== merged cell PLOT ============================= %
%         % Plot fit with data.
%         figure;hold on;
%         xgrid = linspace(0,4,100);
%         histogram(log10(FinalResults{iSamp}.merge_data),xgrid,'Normalization','PDF')
%         h.LineWidth = 2;
%         xlabel( 'PT area (unit: a.u.)', 'Interpreter', 'none' );
%         ylabel( 'PDF', 'Interpreter', 'none' );
    end
end

save_path = 'F:\PROCESS-SPT\2024\20241006-1010_U2OS_TMsignal-OPRM1_PA646_MPALM100';
save(fullfile(save_path,'20241220run_LogDist_2gaussFit.mat'),'data_struct','FinalResults')

%% Plot mobility distribution fitting
close all;
saveDir = 'E:\Submission\2023_MotionBlur\PA646MBX_H2B_1x6xNLS_beads\saved_plots';
% color_str = ["#FF0000", "#0000FF", "#00FF00", "#800080","#FFFF00", "#00FFFF", "#FFA500"];
color_str = lines(length(data_struct));%color_str = prism(2);
titleName = horzcat(data_struct(:).SampleName);
hist_edges = linspace(0,4,100);

if FitParams.DoMergePlot
    % ---------- PDF_pooledFitresults -------------%
    
    titleName = horzcat(data_struct(:).SampleName);
    fig1 = figure('Units','inches','Position',[1.5312 1.8438 4.5 1.5*length(data_struct)]); tiledlayout(length(data_struct),1,'TileSpacing','compact');
    for i  = 1:length(data_struct)
        ax = nexttile(i); hold on;
        
        histogram(log10(FinalResults{i}.merge_data),hist_edges,'Normalization','PDF','FaceColor','none','EdgeColor',color_str(i,:))
        xlabel( 'Log10(PT Area) (unit: a.u.)', 'Interpreter', 'none' );
        ylabel( 'PDF', 'Interpreter', 'none' );
        title(sprintf('%s',titleName(i)),'Interpreter', 'none')
        ax.FontSize = 8;
        ax.XLim = [1.5 4];
    end

    fig2 = figure('Units','inches','Position',[1.5312 1.8438 4.5 1.5*length(data_struct)]); tiledlayout(length(data_struct),1,'TileSpacing','compact');
    
    for i  = 1:length(data_struct)
        ax = nexttile(i); hold on;
        
        histogram(FinalResults{i}.merge_data,0:200:10000,'Normalization','PDF','FaceColor','none','EdgeColor',color_str(i,:))
        xlabel( 'PT area (unit: a.u.)', 'Interpreter', 'none' );
        ylabel( 'PDF', 'Interpreter', 'none' );
        title(sprintf('%s',titleName(i)))
        ax.FontSize = 8;
        % ax.XLim = [0 4000];
        % ax.XLim = [0 200];
    end

    fig33 = figure('Units','inches','Position',[1.5312 1.8438 4.5 1.5*length(data_struct)]); tiledlayout(length(data_struct),1,'TileSpacing','compact');
    for i  = 1:length(data_struct)
        ax = nexttile(i); hold on;

        histogram(log10(FinalResults{i}.merge_D),-3:0.1:2,'Normalization','PDF','FaceColor','none','EdgeColor',color_str(i,:))
        xlabel( 'Log(D)', 'Interpreter', 'none' );
        ylabel( 'PDF', 'Interpreter', 'none' );
        title(sprintf('%s',titleName(i)),'Interpreter', 'none' )
        ax.FontSize = 8;
        ax.XLim = [-3 2];
        % ax.XLim = [-3.5 1];
        % ax.YLim = [0 0.6];
    end

    titleName = horzcat(data_struct(:).SampleName);
    g1 = @(w,x) w(1).*exp(-((x-w(2))./w(3)).^2);
    g2 = @(w,x) w(4).*exp(-((x-w(5))./w(6)).^2);
    g3 = @(w,x) w(1).*exp(-((x-w(2))./w(3)).^2)+w(4).*exp(-((x-w(5))./w(6)).^2);
    fig3 = figure('Units','inches','Position',[1.5312 1.8438 4.5 1.5*length(data_struct)]); tiledlayout(length(data_struct),1,'TileSpacing','compact');
    tilepos = [1 2];
    % fig3 = figure('Units','inches','Position',[1.5312 1.8438 4.5*2 1.5*length(data_struct)/2]); tiledlayout(length(data_struct)/2,2,'TileSpacing','compact');
    for i  = 1:length(data_struct)
        ax = nexttile(i);hold on;
        % ax = nexttile(tilepos(i));hold on;
        xgrid = -3:0.1:2;
        w = coeffvalues(FinalResults{i}.merge_fitResult);
       
            
        histogram(log10(FinalResults{i}.merge_D),-3:0.1:2,'Normalization','PDF','FaceColor','none','EdgeColor',color_str(i,:))
        plot(xgrid,g1(w,xgrid),'--','Color',[0.1 0.1 0.1], 'linewidth', 0.5);
        plot(xgrid,g2(w,xgrid),'--','Color',[0.1 0.1 0.1], 'linewidth', 0.5);
        plot(xgrid,g3(w,xgrid),'k', 'linewidth', 1);
        xlabel( 'Log(D)', 'Interpreter', 'none' );
        ylabel( 'PDF', 'Interpreter', 'none' );
        title(sprintf('%s',titleName(i)))
        ax.FontSize = 8;
        % ax.XLim = [-3.5 1];
        ax.XLim = [-3 2];
        % ax.YLim = [0 0.6];
    end

    fig4 = figure('Units','inches','Position',[1.5312 1.8438 4.5 1.5*length(data_struct)]); tiledlayout(length(data_struct),1,'TileSpacing','compact');
    for i  = 1:length(data_struct)
        ax = nexttile(i); hold on;
        
        histogram(FinalResults{i}.merge_D,0:0.5:20,'Normalization','PDF','FaceColor','none','EdgeColor',color_str(i,:))
        xlabel( 'D', 'Interpreter', 'none' );
        ylabel( 'PDF', 'Interpreter', 'none' );
        title(sprintf('%s',titleName(i)))
        ax.FontSize = 8;
        ax.XLim = [0 20];
    end

    fig5 = figure('Units','inches','Position',[1.5312 1.8438 4.5 1.5]); 
    cmap = jet(length(data_struct));  hold on;
    for i  = 1:length(data_struct)
        % Estimate the PDF using kernel smoothing
        FinalResults{i}.merge_D = FinalResults{i}.merge_D(FinalResults{i}.merge_D ~= 0); % avoid bug
        [f, xi] = ksdensity(log10(FinalResults{i}.merge_D));
        % Plot the kernel smoothed probability
        plot(xi,f,'Color',cmap(i,:),'LineWidth',1)
    end
    xlabel( 'log(D)', 'Interpreter', 'none' );
    ylabel( 'PDF', 'Interpreter', 'none' ); 
    ax = gca; ax.FontSize = 8; ax.XLim = [-3.5 1];% ax.XLim = [-3 2];
    legend(data_struct(:).SampleName,'Location','northwest')
    
    % 3D plot of logD distribution
    data_struct = data_struct([1 2 4 3]);
    fig5 = figure('Units','inches','Position',[1.5312 1.8438 4.5 1.5]); 
    cmap = jet(length(data_struct));  hold on;
    for i  = 1:length(data_struct)
        % Estimate the PDF using kernel smoothing
        FinalResults{i}.merge_D = FinalResults{i}.merge_D(FinalResults{i}.merge_D ~= 0); % avoid bug
        [f, xi] = ksdensity(log10(FinalResults{i}.merge_D));
        % % Plot the kernel smoothed probability
        % plot(xi,f,'Color',cmap(i,:),'LineWidth',1)
    
        % Plot as a line in 3D: X -> class index (k), Y -> log-scaled binCenters, Z -> probabilities
        plot3(xi, repmat(i, size(xi)), f, ...
              'LineWidth', 2, 'Color', cmap(i,:));
    end
    xlabel( 'log(D)', 'Interpreter', 'none' );
    ylabel('Sample')
    zlabel( 'PDF', 'Interpreter', 'none' ); 
    ax = gca; ax.FontSize = 8; ax.XLim = [-3.5 1];% ax.XLim = [-3 2];
    % legend(data_struct(:).SampleName,'Location','northwest')
    view(3)
    grid on
    

    if FitParams.DoMergeFit
        fig5 = figure('Units','inches','Position',[1.5312 1.8438 5.9375 6.2396]); tiledlayout(4,1,'TileSpacing','loose');
        for i  = 1:length(data_struct)
            ax = nexttile(i); hold on;
            histogram(log10(FinalResults{i}.merge_D),-3:0.1:2,'Normalization','PDF','FaceColor','none','EdgeColor',color_str(i,:))
            x_eval = -3:0.1:2;


            y_eval = ft_pdf(FinalResults{i}.w,x_eval);
            plot(x_eval,y_eval,'k', 'linewidth', 0.5);
            xlabel( 'Log10(PT Area) (unit: a.u.)', 'Interpreter', 'none' );
            ylabel( 'PDF', 'Interpreter', 'none' );
            title(sprintf('%s',titleName(i)))
            ax.FontSize = 8;
            ax.XLim = [1.5 4];
%             ax.XLim = [0 200];
        end
        
        
        % heatmap of weights at diffusion coeffcients
        fig6 = figure('Position',[800 420 579 170]);
        temp_merge_norm_w = [];
        for iSamp = 1:length(data_struct)
            temp_merge_norm_w(iSamp,:) = FinalResults{iSamp}.merge_norm_w;
        end
        h = heatmap(temp_merge_norm_w);
        h.ColorLimits = [0 0.15]; % h.ColorScaling = 'scaledrows';
        h.XData = D;
        h.YData = [data_struct(:).SampleName];
        colormap parula
        xlabel('Diffusion coefficient')
        title('Heatmap of normalized weights in each sample')
    end
end


save_path = 'C:\ShareCache\Deng_lab\Shared_denglab_王祖辉\manuscript\MPALM\manuscript_motionBlur\figure_20240328_deleteafteruse';
% exportgraphics(fig1,fullfile(save_path,"LogArea_dist.pdf"))
% exportgraphics(fig2,fullfile(save_path,"Area_dist.pdf"))
exportgraphics(fig3,fullfile(save_path,"RPB1H2B_DMSO_1uMTHZ1_LogD_dist.pdf"))
% exportgraphics(fig4,fullfile(save_path,"D_dist.pdf"))

% % exportgraphics(fig5,fullfile(save_path,"HSF1_HeatShock_LogD_dist.pdf"))
% exportgraphics(fig6,fullfile(save_path,"LogAreaFit_heatmap.pdf"))

%% Locs per frame
% ================== Locs per frame ================== %

x = categorical([data_struct(:).SampleName]);
x = reordercats(x,[data_struct(:).SampleName]);
fig = figure; hold all;
y_mean = cellfun(@(x) mean(x.locPerFrame), FinalResults);
y_std = cellfun(@(x) std(x.locPerFrame), FinalResults);
% y_sem = y_std./sqrt(x_len);
% y_mean = [data_hes_mutant(1:3, 4);merged_model_params(:,3)];
% y_std = [data_hes_mutant(1:3,7);std_model_params(:,3)];

bar(x,y_mean,'FaceColor','flat','EdgeColor','none');  % categorical(label) is the x label
ylabel('# locs per frame per nuclues', 'FontSize', 20); 
% overlay with error bar
errorbar(x,y_mean, y_std, '.', 'color', 'k', 'linewidth', 1); hold off % errorbar(x,y,err)

% overlay data point
hold on
xCenter = 1:length(FinalResults); 
spread = 0.5; % 0=no spread; 0.5=random spread within box bounds (can be any value)
for i = 1:length(FinalResults)
    allData = FinalResults{i}.locPerFrame; % bound fraction
    plot(rand(size(allData))*spread -(spread/2) + xCenter(i), allData, 'k.','MarkerSize', 5)
end

% adjust format
ax = gca; ax.FontSize = 12; ax.LineWidth = 0.5; ax.XAxis.TickLabelRotation = 45;  ax.TickLabelInterpreter = 'none';  ax.XColor = 'k'; ax.YColor = 'k';   %Relative length of each axis, specified as a three-element vector of the form [px py pz] 
ax.PlotBoxAspectRatio = [1 1.5 1]; %ax.YLim = [0, 0.8]; 
exportgraphics(fig,fullfile("dense_sparse_locDensity_PerNucleus.pdf"))



%% Run RDF
%close all;clc;
draw_roi = false;
load_roi = false;
keep_same_density = false;
filt_D = false;

minRawLocDensity = 5.4;

% loop through each sample
for iSamp = 3%1:length(data_struct)
    fprintf('Loading sample: %s \n',data_struct(iSamp).SampleName);
    tab_all_SR_pred = table();
    temp_single_data = {};
    tot_cell_num = length(data_struct(iSamp).workspaces);    
    data_struct(iSamp).MeanGrScore = [];
    data_struct(iSamp).MeanGr0Score = [];
    data_struct(iSamp).GrScore = [];
    data_struct(iSamp).Gr0Score = [];

    for iCell = 1:tot_cell_num
        fprintf('Process cell: %s \n',data_struct(iSamp).workspaces(iCell).name);
        clear tab_SR_pred roiVertices params;       
        %close all;

        % RDF calculation params 
        params.dr = 10;
        params.calib = 30;
        params.dxy = 30;
        params.camPixelSize_Ch1 = 110;
        params.num_pts = 30;
        params.D_SR_const = 2007;
        params.maxD = 0.05;
        
        SRFile = dir(fullfile(data_struct(iSamp).workspaces(iCell).folder, data_struct(iSamp).workspaces(iCell).name, 'UNet_*.csv'));
        tab_SR_pred = readtable(fullfile(SRFile.folder, SRFile.name));

        if filt_D
           tab_SR_pred.D = (tab_SR_pred.SRArea./params.D_SR_const).^2;
           tab_SR_pred = tab_SR_pred(tab_SR_pred.D<=params.maxD,:);
        end

        if draw_roi
            fig = figure('Position', [471 104 822 738]);
            
            scatter(tab_SR_pred.Xpos, tab_SR_pred.Ypos, 2, 'filled', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.3);
            xlabel('X Position');
            ylabel('Y Position');
            title('Select ROI');
            axis image ij
    
            % Draw a polygon ROI
            hROI = drawpolygon('LineWidth', 2, 'Color', 'r');
    
            % Get the X and Y coordinates of the vertices of the polygon
            roiVertices = hROI.Position;
    
            % Use inpolygon to find points within the ROI
            [in, on] = inpolygon(tab_SR_pred.Xpos, tab_SR_pred.Ypos, roiVertices(:, 1), roiVertices(:, 2));
    
            % Filter the table to include only rows within the ROI
            tab_SR_pred = tab_SR_pred(in | on, :);
            close(fig)
        end
        
        % tab_SR_pred.Xpos/Ypos unit is camera pixel
        if load_roi

            load(fullfile(SRFile.folder,sprintf('roi_metadata.mat')),'hROI');

            % Get the X and Y coordinates of the vertices of the polygon
            roiVertices = hROI.Position;

            fig = figure('Position', [471 104 822 738]);            
            scatter(tab_SR_pred.Xpos, tab_SR_pred.Ypos, 2, 'filled', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.3);
            xlabel('X Position');
            ylabel('Y Position');
            title('Select ROI');
            axis image ij
    
            % Draw a polygon ROI
            hROI = drawpolygon('Position',roiVertices,'LineWidth', 2, 'Color', 'r');
                    
            % Use inpolygon to find points within the ROI
            [in, on] = inpolygon(tab_SR_pred.Xpos, tab_SR_pred.Ypos, roiVertices(:, 1), roiVertices(:, 2));
    
            % Filter the table to include only rows within the ROI
            tab_SR_pred = tab_SR_pred(in | on, :);

            
            close(fig)
        end

        if ~draw_roi && ~load_roi
            % Assume tab_SR_pred is a table with Xpos and Ypos columns
            x = tab_SR_pred.Xpos;
            y = tab_SR_pred.Ypos;
            
            % Calculate the minimum and maximum coordinates
            minX = min(x);
            maxX = max(x);
            minY = min(y);
            maxY = max(y);
            
            % Define the rectangle ROI covering all points
            roiVertices = [minX, minY; minX , maxY; maxX , maxY;maxX, minY];           
        end

        if keep_same_density
                       
            % Calculate the area of the polygon ROI using the polyarea function
            roiArea = polyarea(roiVertices(:, 1), roiVertices(:, 2));           

            % Calculate the number of rows to select
            numRowsToSelect = round(minRawLocDensity*roiArea);
            
            % Generate a random permutation of row indices
            randomIndices = randperm(height(tab_SR_pred), numRowsToSelect);
            
            % Select the rows using the generated indices
            tab_SR_pred = tab_SR_pred(randomIndices, :);

        end

        % change Xpos/Ypos unit to nm
        roiVertices = roiVertices*params.camPixelSize_Ch1;
        tab_SR_pred.Xpos = tab_SR_pred.Xpos*params.camPixelSize_Ch1;
        tab_SR_pred.Ypos = tab_SR_pred.Ypos*params.camPixelSize_Ch1;

        Blurdata_mat = dir(fullfile(data_struct(iSamp).workspaces(iCell).folder, data_struct(iSamp).workspaces(iCell).name, 'Blurdata_*.mat'));
        Blurdata_mat = load(fullfile(Blurdata_mat(1).folder, Blurdata_mat(1).name), 'impars');
    
        params.ImWidth = Blurdata_mat.impars.ImWidth;
        params.ImHeight = Blurdata_mat.impars.ImHeight;
                          
        % Calculate the mean g(r) score for the filtered localizations
        % RDF calculation params        
        params.Ch1_loc_XYLim = [0,params.camPixelSize_Ch1*params.ImWidth,0,params.camPixelSize_Ch1*params.ImHeight];
        [MeanGrScore,MeanGr0Score,GrScore,Gr0Score,params] = GetROIMeanGrScore(tab_SR_pred, roiVertices, params);
        data_struct(iSamp).MeanGrScore = [data_struct(iSamp).MeanGrScore;MeanGrScore];
        data_struct(iSamp).MeanGr0Score = [data_struct(iSamp).MeanGr0Score;MeanGr0Score];
        data_struct(iSamp).GrScore = [data_struct(iSamp).GrScore;GrScore];
        data_struct(iSamp).Gr0Score = [data_struct(iSamp).Gr0Score;Gr0Score];
    end    
end

% save(fullfile(rootPath,'RDF_data_struct.mat'),'data_struct','params')


% Plot the result

figure('Units','inches','Position',[7.1562 3.8854 5.8333 2.7812]);
color_str = lines(4);%color_str = prism(2);
titleName = horzcat(data_struct(:).SampleName);
hold on;
for i = 3%1:length(data_struct)   
    errorbar(params.radii,mean(data_struct(i).MeanGrScore,1,"omitnan"),std(data_struct(i).MeanGrScore,0,1,"omitnan"),'s-','MarkerFaceColor',color_str(i,:),'MarkerEdgeColor',color_str(i,:),'Color',color_str(i,:),'MarkerSize',2);
end
plot([0 max(params.radii)], [1 1], 'k--', 'LineWidth', 1, 'HandleVisibility','off');
ylabel('G(r)');xlabel('r(nm)');
% ylim([0 12])
% xlim([0 300])
legend(titleName)



%% Run RDP version2
close all;clc;
% RDF calculation params 
% params.dr = 30;
% params.calib = 0.5*params.dr;
params.dxy = 30; % pixel size of Im2 used to calculate auto-correlation, unit nm
params.camPixelSize_Ch1 = 110;
params.maxD = 0.05;
params.D_SR_const = 1162;
% params.num_pts = 15;

draw_roi = true;
load_roi = false;
keep_same_density = false;minRawLocDensity = 0.1;
filt_D = false;
filt_frame = false;


% if keep_same_density
%     for iSamp = 1:length(data_struct)
%     fprintf('Loading sample: %s \n',data_struct(iSamp).SampleName);
%     tab_all_SR_pred = table();
%     temp_single_data = {};
%     tot_cell_num = length(data_struct(iSamp).workspaces);    
%     data_struct(iSamp).MeanGrScore = [];
%     data_struct(iSamp).MeanGr0Score = [];
%     data_struct(iSamp).RawLocDensity = [];
% 
%         for iCell = 1:tot_cell_num
%             fprintf('Process cell: %s \n',data_struct(iSamp).workspaces(iCell).name);
%             clear tab_SR_pred;
%             close all;
%             SRFile = dir(fullfile(data_struct(iSamp).workspaces(iCell).folder, data_struct(iSamp).workspaces(iCell).name, 'UNet_*.csv'));
%             tab_SR_pred = readtable(fullfile(SRFile.folder, SRFile.name));
% 
%             if load_roi
%                 load(fullfile(SRFile.folder,sprintf('roi_metadata.mat')),'hROI');
% 
%                 % Get the X and Y coordinates of the vertices of the polygon
%                 roiVertices = hROI.Position;
% 
%                 fig = figure('Position', [471 104 822 738]);            
%                 scatter(tab_SR_pred.Xpos, tab_SR_pred.Ypos, 2, 'filled', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.3);
%                 xlabel('X Position');
%                 ylabel('Y Position');
%                 title('Select ROI');
%                 axis image ij
% 
%                 % Draw a polygon ROI
%                 hROI = drawpolygon('Position',roiVertices,'LineWidth', 2, 'Color', 'r');               
% 
%                 % Use inpolygon to find points within the ROI
%                 [in, on] = inpolygon(tab_SR_pred.Xpos, tab_SR_pred.Ypos, roiVertices(:, 1), roiVertices(:, 2));
% 
%                 % Filter the table to include only rows within the ROI
%                 tab_SR_pred = tab_SR_pred(in | on, :);                                  
% 
%                 % Calculate the number of points inside the ROI
%                 numPointsInROI = sum(in);
% 
%                 % Calculate the area of the polygon ROI using the polyarea function
%                 roiArea = polyarea(roiVertices(:, 1), roiVertices(:, 2));
% 
%                 % Compute the localization density (number per unit area)
%                 localizationDensity = numPointsInROI / roiArea;                
% 
%                 close(fig)
%             end
% 
%             data_struct(iSamp).RawLocDensity(iCell) = localizationDensity;
%         end
%     end
% end
% 
% minRawLocDensity = cellfun(@min,{data_struct(:).RawLocDensity});
% minRawLocDensity = min(minRawLocDensity);




% loop through each sample
for iSamp = 1:length(data_struct)
    fprintf('Loading sample: %s \n',data_struct(iSamp).SampleName);
    tab_all_SR_pred = table();
    temp_single_data = {};
    tot_cell_num = length(data_struct(iSamp).workspaces);    
    data_struct(iSamp).MeanGrScore = [];
    data_struct(iSamp).MeanGr0Score = [];
    data_struct(iSamp).GrScore = [];

    for iCell = 1:tot_cell_num
        fprintf('Process cell: %s \n',data_struct(iSamp).workspaces(iCell).name);
        clear tab_SR_pred;
        close all;
        SRFile = dir(fullfile(data_struct(iSamp).workspaces(iCell).folder, data_struct(iSamp).workspaces(iCell).name, 'drift_*.csv'));
        tab_SR_pred = readtable(fullfile(SRFile.folder, SRFile.name));

        if filt_frame           
           tab_SR_pred = tab_SR_pred(tab_SR_pred.Frame<=1000,:);
        end

        if filt_D
           tab_SR_pred.D = (tab_SR_pred.SRArea./params.D_SR_const).^2;
           tab_SR_pred = tab_SR_pred(tab_SR_pred.D<=params.maxD,:);
        end

        if draw_roi
            fig = figure('Position', [471 104 822 738]);
            
            scatter(tab_SR_pred.Xpos, tab_SR_pred.Ypos, 2, 'filled', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.3);
            xlabel('X Position');
            ylabel('Y Position');
            title('Select ROI');
            axis image ij
    
            % Draw a polygon ROI
            hROI = drawpolygon('LineWidth', 2, 'Color', 'r');
    
            % Get the X and Y coordinates of the vertices of the polygon
            roiVertices = hROI.Position;
    
            % Use inpolygon to find points within the ROI
            [in, on] = inpolygon(tab_SR_pred.Xpos, tab_SR_pred.Ypos, roiVertices(:, 1), roiVertices(:, 2));
    
            % Filter the table to include only rows within the ROI
            tab_SR_pred = tab_SR_pred(in | on, :);
            close(fig)
        end
        
        % tab_SR_pred.Xpos/Ypos unit is camera pixel
        if load_roi

            load(fullfile(SRFile.folder,sprintf('roi_metadata.mat')),'hROI');

            % Get the X and Y coordinates of the vertices of the polygon
            roiVertices = hROI.Position;

            fig = figure('Position', [471 104 822 738]);            
            scatter(tab_SR_pred.Xpos, tab_SR_pred.Ypos, 2, 'filled', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.3);
            xlabel('X Position');
            ylabel('Y Position');
            title('Select ROI');
            axis image ij
    
            % Draw a polygon ROI
            hROI = drawpolygon('Position',roiVertices,'LineWidth', 2, 'Color', 'r');
                    
            % Use inpolygon to find points within the ROI
            [in, on] = inpolygon(tab_SR_pred.Xpos, tab_SR_pred.Ypos, roiVertices(:, 1), roiVertices(:, 2));
    
            % Filter the table to include only rows within the ROI
            tab_SR_pred = tab_SR_pred(in | on, :);
            
            close(fig)
        end

        if keep_same_density
                       
            % Calculate the area of the polygon ROI using the polyarea function
            roiArea = polyarea(roiVertices(:, 1), roiVertices(:, 2));           

            % Calculate the number of rows to select
            numRowsToSelect = round(minRawLocDensity*roiArea);
            
            % Generate a random permutation of row indices
            randomIndices = randperm(height(tab_SR_pred), numRowsToSelect);
            
            % Select the rows using the generated indices
            tab_SR_pred = tab_SR_pred(randomIndices, :);

        end

        
        roiVertices = roiVertices*params.camPixelSize_Ch1;
        

        Blurdata_mat = dir(fullfile(data_struct(iSamp).workspaces(iCell).folder, data_struct(iSamp).workspaces(iCell).name, 'Blurdata_*.mat'));
        Blurdata_mat = load(fullfile(Blurdata_mat(1).folder, 'Blurdata_UNet_mask_MBX_20240620_2035_epoch20_Ch1.mat'), 'impars');
    
        params.ImWidth = Blurdata_mat.impars.ImWidth;
        params.ImHeight = Blurdata_mat.impars.ImHeight;

        tab_SR_pred.Xpos = tab_SR_pred.Xpos*params.camPixelSize_Ch1;
        tab_SR_pred.Ypos = tab_SR_pred.Ypos*params.camPixelSize_Ch1;

        % generate density based PALM image                   
        loc_XYLim = [0,params.camPixelSize_Ch1*params.ImWidth,0,params.camPixelSize_Ch1*params.ImHeight];
        Xmin = loc_XYLim(1);Xmax=loc_XYLim(2);Ymin = loc_XYLim(3);Ymax=loc_XYLim(4);                            
        Edges{1}=Xmin:params.dxy:Xmax+params.dxy; % make edge cover all data points
        Edges{2}=Ymin:params.dxy:Ymax+params.dxy;                    
        [Im,~,~,~,~] = histcounts2(tab_SR_pred.Xpos,tab_SR_pred.Ypos,Edges{1},Edges{2}); % bin every localizations into grid, similar to Im = hist3([Xpos',Ypos'],'Edges',Edges);                   
        Im2 = transpose(Im);     
        Im3 = Im2 > 1;
        Im3 = double(Im3);

        
            
        xy_position = [tab_SR_pred.Xpos, tab_SR_pred.Ypos];
    
        % Set the boundaries for random rectangle generation
        max_x = ceil(max(xy_position(:, 1)));
        max_y = ceil(max(xy_position(:, 2)));
        min_x = floor(min(xy_position(:, 1)));
        min_y = floor(min(xy_position(:, 2)));
        total_mask_num = 1;
        mask_num = total_mask_num;
        g_r = [];              

        while mask_num > 0    
            % Generate random coordinates within the ROI bounds
            x = randi([min_x+500, max_x-500]); % Adjusted to prevent exceeding bounds, unit nm
            y = randi([min_y+500, max_y-500]);  
    
            % Define a rectangle (adjust size as needed), unit nm
            w = 4000;
            h = 4000;
    
            % Check if the entire rectangle is within the ROI
            if all(inpolygon([x, x, x+w, x+w], [y, y+h, y, y+h], roiVertices(:, 1), roiVertices(:, 2)))
                % Collect localizations within the random rectangle
                idx = xy_position(:, 1) >= x & xy_position(:, 1) <= x + w & ...
                      xy_position(:, 2) >= y & xy_position(:, 2) <= y + h;
                locs_in_rect = xy_position(idx, :);
    
                % Calculate properties for g(r) using locs_in_rect
                if ~isempty(locs_in_rect)                                                                                       
                    % figure;
                    % imshow(Im2,[0 1])
                    % ax = gca;
                    % % hROI = drawrectangle('LineWidth', 2, 'Color', 'r'); % Interactively draw an ROI
                    % title('Select mask');
                    % 
                    % 
                    % if mask_num == total_mask_num
                    %     previousROI = [floor(x/params.dxy),floor(y/params.dxy),w/params.dxy,h/params.dxy];                                          
                    % end
                    % hROI = drawrectangle(ax,'Position',previousROI,'LabelVisible','on',...
                    %                 'FaceAlpha', 0,'LineWidth', 2,'Color','green',...
                    %                 'InteractionsAllowed','translate');
                    % f = msgbox("Finish ROI adjustment?"); f.Position = [921.7500 406.5000 150 51.7500];
                    % uiwait(f);
                    % previousROI = hROI.Position;
                    % mask = createMask(hROI);


                    % Initialize the mask as a logical matrix
                    mask = false(size(Im3));
                    % Set the rectangle area to true
                    mask(floor(y/params.dxy):ceil((y+h)/params.dxy), floor(x/params.dxy):ceil((x+w)/params.dxy)) = true;       

                    % Overlay ROI on rendered image (debug only)
                    figure;
                    imshow(Im2,[0 1]);hold on;
                    B = bwboundaries(mask);
                    boundary = B{1};
                    plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2); hold off;
                                       
                    fprintf("Calculating RDF for the selected ROI ...count down %d \n",mask_num);
                                                    
                    % Get the mask from the drawn ROI                                
                    [G, r, g, dg, mask] = get_autocorr(Im3 , mask, 450, true);
                    g_r = [g_r;g(2:length(r))];

                    % % Fitting to equation 2 (random model)
                    % sigma_s = 1;
                    % rho_p = 5.4;
                    % gr_eq2 = @(r,ss,rho_p) 1/(4*pi*ss*rho_p)*exp(-r^2/(4*ss))+1;
                    % 
                    % y = g(2:length(r));
                    % x = r(2:length(r));
                    % 
                    % % ---------- Non-random fitting model ----------
                    % % Initial guess for parameter s
                    % initialS = [1,1,1]; % theta = param(1); ss= param(2); rho_p= param(3);
                    % 
                    % % Fit the model to the data
                    % [sFit, resnorm] = lsqcurvefit(@customModel, initialS, x, y);
                    % 
                    % % Display the fitted parameter
                    % disp(['Fitted parameter s: ', num2str(sFit)]);
                    % 
                    % % Generate fitted curve
                    % yFitted = customModel(sFit, x);
                    % 
                    % % Plot original data and fitted curve
                    % figure;
                    % plot(x, y, 'bo', 'DisplayName', 'Original Data');
                    % hold on;
                    % plot(x, yFitted, 'r-', 'DisplayName', 'Fitted Curve');
                    % xlabel('x');
                    % ylabel('y');
                    % legend;
                    % title('Data Fitting with Custom Model');

                    mask_num = mask_num - 1;
                end
            end
        end
        MeanGrScore = mean(g_r,1);
        data_struct(iSamp).MeanGrScore = [data_struct(iSamp).MeanGrScore;MeanGrScore];
        data_struct(iSamp).GrScore = [data_struct(iSamp).GrScore;g_r];
    end    
end
params.radii = r(2:length(r))*params.dxy;


% Plot the result

figure('Units','inches','Position',[7.1562 3.8854 5.8333 2.7812]);
color_str = lines(4);%color_str = prism(2);
titleName = horzcat(data_struct(:).SampleName);
hold on;
for i = 1:length(data_struct)   
    % plot(params.radii,data_struct(i).MeanGrScore,'Color',color_str(i,:))
    % errorbar(params.radii,mean(data_struct(i).MeanGrScore,1,"omitnan"),std(data_struct(i).MeanGrScore,0,1,"omitnan"),'s-','MarkerFaceColor',color_str(i,:),'MarkerEdgeColor',color_str(i,:),'Color',color_str(i,:),'MarkerSize',2);
    
    % Errorbar use SEM
    % errorbar(params.radii,mean(data_struct(i).MeanGrScore,1,"omitnan"),std(data_struct(i).MeanGrScore,0,1,"omitnan")/sqrt(height(data_struct(i).MeanGrScore)),'s-','MarkerFaceColor',color_str(i,:),'MarkerEdgeColor',color_str(i,:),'Color',color_str(i,:),'MarkerSize',2);
    shadedErrorBar(params.radii,mean(data_struct(i).MeanGrScore,1,"omitnan"),std(data_struct(i).MeanGrScore,0,1,"omitnan")/sqrt(height(data_struct(i).MeanGrScore)),'lineProps',{'-or','color',color_str(i,:)})
                                
end
plot([0 max(params.radii)], [1 1], 'k--', 'LineWidth', 1, 'HandleVisibility','off');
ylabel('G(r)');xlabel('r(nm)');
% ylim([0 12])
xlim([0 100])
legend(titleName)
axis square

%% test fitting Sengupta 2013 Nat Met.

clc;

% ---------- Non-random fitting model ----------
% Initial guess for parameter s
initialS = [1000,10,1,1]; % A = param(1); theta = param(2); ss= param(3); rho_p= param(4);

% Fit the model to the data
[sFit, resnorm] = lsqcurvefit(@customModel, initialS, x, y);

% Generate fitted curve
yFitted = customModel(sFit, x);

% Plot original data and fitted curve
figure;
plot(x, y, 'bo', 'DisplayName', 'Original Data');
hold on;
plot(x, yFitted, 'r-', 'DisplayName', 'Fitted Curve');
xlabel('x');
ylabel('y');
legend;
title('Data Fitting with Custom Model');




%% simulation RDF test

x_min = 10;   % Minimum x-coordinate (left edge)
y_min = 20;   % Minimum y-coordinate (bottom edge)
width = 100;   % Width of the rectangle
height = 100;  % Height of the rectangle

% Number of points to generate
num_points = 10000;

% Generate random x and y coordinates within the rectangle
x_random = x_min + width * rand(num_points, 1);
y_random = y_min + height * rand(num_points, 1);

% Define parameters
params.calib = 1;
params.dr = 0.1;
params.dxy = params.dr;

% params.camPixelSize_Ch1 = 110;
params.num_pts = 15;

% % Example random points (replace with your actual data)
% x_random = rand(1000, 1) * 50 + 10; % Random X coordinates
% y_random = rand(1000, 1) * 30 + 20; % Random Y coordinates

% Define grid edges based on your coordinate ranges
x_edges = min(x_random):params.dxy:max(x_random)+params.dxy;
y_edges = min(y_random):params.dxy:max(y_random)+params.dxy;

data_x = (x_random-min(x_random))./params.dxy+0.5;
data_y = (y_random-min(y_random))./params.dxy+0.5;

% Bin localizations into a grid
[N, ~, ~,~,~] = histcounts2(x_random, y_random, x_edges, y_edges);

% Display the binned image
figure; 
imshow(N,[],'InitialMagnification',1000);hold on;
% imagesc(x_edges, y_edges, N);
scatter(data_y,data_x,8, 'filled')
title('Binned Localizations');
hROI = drawrectangle('LineWidth', 2, 'Color', 'r'); % Interactively draw an ROI

% Get the mask from the drawn ROI
mask = createMask(hROI); % Logical mask of the selected ROI

fprintf("Calculating RDF for the selected ROI ...\n");

% % -------------------- Option 1
% % Calculate autocorrelation using binned data inside the ROI
% [corrs, params] = ParCalculateXcorrApp(data_x, data_y, mask, params);
% 
% figure('Units','inches','Position',[7.1562 3.8854 5.8333 2.7812]);
% hold on;
% errorbar(params.radii,corrs.C_11,corrs.dC_11,'s-','MarkerFaceColor','b','MarkerEdgeColor','b','Color','b','MarkerSize',2)
% plot([0 max(params.radii)], [1 1], 'k--', 'LineWidth', 1, 'HandleVisibility','off');
% ylabel('G(r)');xlabel('r(nm)');
% ylim([0 12])
% %xlim([0 300])
% axis square



% -------------------- Option 2
[G, r, g, dg, mask] = get_autocorr(N , mask, 50,true);




%% (legacy) U-Net area 
% loop through each sample
clear FinalResults FinalResults
for iSamp = 1:length(data_struct)

    % loop through each cell   
    tab_all_SR_pred = table();
    tot_cell_num = length(data_struct(iSamp).workspaces);
    for iCell = 1:tot_cell_num
        clear tab_SR_pred
        tab_SR_pred = readtable(fullfile(data_struct(iSamp).workspaces(iCell).folder,data_struct(iSamp).workspaces(iCell).name, data_struct(iSamp).SRFile));
        
        fprintf('localization number %d ...\n',height(tab_SR_pred))
        tab_all_SR_pred = [tab_all_SR_pred;tab_SR_pred];
        
    end
    FinalResults{iSamp}.merge_data = tab_all_SR_pred.UNetArea;
%         clear tab_psf_fitresult
%         load(fullfile(data_struct(iSamp).workspaces(iCell).folder,data_struct(iSamp).workspaces(iCell).name,data_struct(iSamp).UNetMat),'tab_psf_fitresult');
%         fprintf('localization number %d ...\n',height(tab_psf_fitresult))
%         tab_all_psf_result = [tab_all_psf_result;tab_psf_fitresult];
%         FinalResults{iSamp}.merge_data = tab_all_psf_result.Area;
    
end

saveDir = 'E:\Submission\2023_MotionBlur\PA646MBX_H2B_1x6xNLS_beads\saved_plots';
% color_str = ["#0072BD","#D95319","#7E2F8E","#77AC30","#F9F1A5","#000000"];
color_str = ["#FF0000", "#0000FF", "#00FF00", "#FFFF00", "#800080", "#00FFFF", "#FFA500"];

if FitParams.DoMergePlot
    % ---------- PDF_pooledFitresults -------------%
    
    titleName = horzcat(data_struct(:).SampleName);   
    fig = figure('Units','inches','Position',[1.5312 1.8438 5.9375 1.5*length(data_struct)]); tiledlayout(length(data_struct),1,'TileSpacing','compact');
    
    for i  = 1:length(data_struct)
        ax = nexttile(i); hold on;
        histogram(FinalResults{i}.merge_data,20:5:200,'Normalization','PDF','FaceColor','none','EdgeColor',color_str(i,:))
        xlabel( 'UNet area (unit: pix)', 'Interpreter', 'none' );
        ylabel( 'PDF', 'Interpreter', 'none' );
        title(sprintf('%s',titleName(i)))
        ax.FontSize = 8;
        ax.XLim = [0 200];
    end
end
    
%% Auxiliary functions

function [x,resnorm,residual,exitflag,output] = multi_gauss_fit(data,ft)
xgrid = linspace(0,4,100);
pd = fitdist(data,'Kernel');
y = cdf(pd,xgrid);

[xData, yData] = prepareCurveData( xgrid, y );

% Fit model to data.
x0 = [repmat(0.0556,[1 18]), ones(1,18)*0.2]; % Initial values for w and sigma
lb = [zeros(1,18), zeros(1,18)]; % Lower bounds for w and sigma
ub = [Inf(1,18), Inf(1,18)]; % Upper bounds for w and sigma

% sigma_sim = [0.264707351793163	0.216854887410073	0.237400174745841	0.219043205025058	0.219676454766064	0.226765723973284	0.229884139969482	0.235904013254032	0.240858161844799	0.218397291217037	0.195818439267191	0.199787464542736	0.196803573519027	0.183287247164872	0.168473705433979	0.135959323108265	0.100974382396461	0.110644103340083];
% x0 = [repmat(0.0556,[1 18]), sigma_sim]; % Initial values for w and sigma
% lb = [zeros(1,18), sigma_sim.*0.1]; % Lower bounds for w and sigma
% ub = [Inf(1,18), sigma_sim.*10]; % Upper bounds for w and sigma

opts = optimset('MaxIter',5000,'MaxFunEvals', 5000, 'TolFun',1e-6,'TolX',1e-5); % best performance for benchmark molecules
[x,resnorm,residual,exitflag,output] = lsqcurvefit(ft,x0,xData,yData,lb,ub,opts);
fprintf('finish fitting\n')
end

function [x,resnorm,residual,exitflag,output] = multi_gauss_fit_fixSigma(data,ft)
xgrid = linspace(0,4,100);
pd = fitdist(data,'Kernel');
y = cdf(pd,xgrid);

[xData, yData] = prepareCurveData( xgrid, y );

% opts = fitoptions('Method','NonlinearLeastSquares', ...
%     'Lower',zeros(1,12), ...
%     'Upper',ones(1,12), ...
%     'StartPoint',repmat(0.0833,[1 12]), ... 0.045 equals 1/12 do not use others
%     'DiffMaxChange',1e-3,...
%     'TolFun',1e-2,... default value 1e-6 gives more smoothed weights but Fitting not always converged to a solution
%     'Display','final');

% Fit model to data.
x0 = repmat(0.0555,[1 18]);
lb = zeros(1,18);
ub = Inf(1,18);
opts = optimset('MaxIter',5000,'MaxFunEvals', 5000, 'TolFun',1e-20,'TolX',1e-20);

[x,resnorm,residual,exitflag,output] = lsqcurvefit(ft,x0,xData,yData,lb,ub,opts);
fprintf('finish fitting\n')
end

function fitresult = cell_fit_gaussian2(xData,yData,zData,pixel_size)
%cell_fit_gaussian2 Fitting beads PSF with 2D symmetric Gaussian
%   Detailed explanation goes here

% assign start point
EmissionWavelength = 664;
impars.PixelSize=pixel_size; % um per pixel, Andor iXon Ultra 897
impars.psf_scale=1.35; % PSF scaling
impars.wvlnth= EmissionWavelength/1000; %emission wavelength in um
impars.NA=1.49; % NA of detection objective
impars.psfStd= impars.psf_scale*0.55*(impars.wvlnth)/impars.NA/1.17/impars.PixelSize/2; % PSF standard deviation in pixels

% assign the fit model using MTT paper
g = fittype('I * exp(-(1/(2*r0^2))*((y - y0)^2 + (x - x0)^2)) + m','independent',{'x','y'},'dependent','z');
fo = fitoptions(g);
fo.StartPoint = [max(zData),min(zData),impars.psfStd+1,mean(xData),mean(yData)];
fo.Lower = [0,0,1.40,0,0]; 
fo.Upper = [Inf,Inf,Inf,32,32];
fo.Robust = 'Bisquare';
fo.MaxIter = 1000;
fo.TolFun = 1.0000e-08;
fo
[fitobject,gof]  = fit([xData,yData],zData,g,fo);
gof
fitresult = [fitobject.I, fitobject.m,fitobject.r0, fitobject.x0,fitobject.y0];
end

function [fitresult, gof,output] = two_gauss_fit(data)
% Compute the kernel density estimate

% For protein in nucleus
xgrid = -3:0.1:2; 

% % For protein on cytoplasma membrane
% xgrid = -3.5:0.1:1; 


% xgrid = linspace(min(data),max(data),100);
pd = fitdist(data,'Kernel');
y = pdf(pd,xgrid);

[xData, yData] = prepareCurveData( xgrid, y );

% Set up fittype and options.
ft = fittype( 'gauss2' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares','Algorithm','Trust-Region');
opts.Display = 'Off';

% For protein in nucleus
opts.Lower = [0 -3 0 0 -0.3010 0];
opts.StartPoint = [0.5 -1.6 0.1 0.5 0.6 0.1]; 
opts.Upper = [Inf -1.3010 1 Inf 1.47 1];

% % For protein on cytoplasma membrane
% opts.Lower = [0 -2.5 0 0 -1 0];
% opts.StartPoint = [0.5 -2 0.1 0.5 0 0.1]; 
% opts.Upper = [Inf -1.5 0.7 Inf 0.5 1];

% Fit model to data.
[fitresult, gof, output] = fit( xData, yData, ft, opts );
end


function [MeanGrScore,MeanGr0Score,GrScore,Gr0Score,params] = GetROIMeanGrScore(tab_SR_pred_filtered, roiVertices, params)
    % GetROIMeanGrScore Summary of this function goes here
    % Detailed explanation goes here

    ImOut_density = density_PALM(tab_SR_pred_filtered,tab_SR_pred_filtered.Xpos,tab_SR_pred_filtered.Ypos,'all',params.Ch1_loc_XYLim);

    xy_position = [tab_SR_pred_filtered.Xpos, tab_SR_pred_filtered.Ypos];
    
    % Set the boundaries for random rectangle generation
    max_x = ceil(max(xy_position(:, 1)));
    max_y = ceil(max(xy_position(:, 2)));
    min_x = floor(min(xy_position(:, 1)));
    min_y = floor(min(xy_position(:, 2)));
    
    total_mask_num = 3;
    mask_num = total_mask_num;
    g_r = [];
    g_r0 = []; % g_r for uniform background with same density as sample roi    

    while mask_num > 0    
        % Generate random coordinates within the ROI bounds
        x = randi([min_x+1000, max_x-1000]); % Adjusted to prevent exceeding bounds, unit nm
        y = randi([min_y+1000, max_y-1000]);  

        % Define a rectangle (adjust size as needed), unit nm
        w = 4000;
        h = 4000;

        % Check if the entire rectangle is within the ROI
        if all(inpolygon([x, x, x+w, x+w], [y, y+h, y, y+h], roiVertices(:, 1), roiVertices(:, 2)))
            % Collect localizations within the random rectangle
            idx = xy_position(:, 1) >= x & xy_position(:, 1) <= x + w & ...
                  xy_position(:, 2) >= y & xy_position(:, 2) <= y + h;
            locs_in_rect = xy_position(idx, :);

            % Calculate properties for g(r) using locs_in_rect
            if ~isempty(locs_in_rect)                

                Xmin = params.Ch1_loc_XYLim(1);Xmax=params.Ch1_loc_XYLim(2);Ymin = params.Ch1_loc_XYLim(3);Ymax=params.Ch1_loc_XYLim(4);
                data1_x = (tab_SR_pred_filtered.Xpos-Xmin)./params.dxy+0.5 ; % unit rendered pixel
                data1_y = (tab_SR_pred_filtered.Ypos-Ymin)./params.dxy+0.5 ; % unit rendered pixel

                
                % Overlay ROI on rendered image (debug only)
                figure;
                imshow(ImOut_density,[0 5000],'InitialMagnification',1000,'colormap',hot);ax = gca;
                % hold on;
                % B = bwboundaries(mask);
                % boundary = B{1};
                % plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2); hold off;

                if mask_num == total_mask_num
                    previousROI = [floor(x/params.dxy),floor(y/params.dxy),w/params.dxy,h/params.dxy];                                          
                end
                hROI = drawrectangle(ax,'Position',previousROI,'LabelVisible','on',...
                                'FaceAlpha', 0,'LineWidth', 2,'Color','green'); %'InteractionsAllowed','translate'
                f = msgbox("Finish ROI adjustment?"); f.Position = [921.7500 406.5000 150 51.7500];
                uiwait(f);
                previousROI = hROI.Position;
                mask = createMask(hROI);
    
                % % Initialize the mask as a logical matrix
                % mask = false(size(ImOut_density));                
                % % Set the rectangle area to true
                % mask(floor(y/params.dxy):ceil((y+h)/params.dxy), floor(x/params.dxy):ceil((x+w)/params.dxy)) = true;
                
                % roi = drawrectangle(app.densityPALMfig.Children(2));  
                % [m,n] = size(ImOut_density);
                % bw = createMask(roi,m,n);                
                
                fprintf("Calculating RDF for the selected ROI ...count down %d \n",mask_num);
                [corrs, params] = ParCalculateXcorrApp(data1_x, data1_y, mask, params);                 
                

                

                % Generate simulated image with same but uniform density
                % distribution to calculate g_r

                % estimate loc density in roi, unit loc/rendered_pixel^2
                ave_density = height(locs_in_rect)/sum(sum(mask));

                x_min = 0;   % Minimum x-coordinate (left edge)
                y_min = 0;   % Minimum y-coordinate (bottom edge)
                Width = 500;  % Width of the rectangle
                Height = 500; % Height of the rectangle
                
                % Number of points to generate
                num_points = ave_density*Width*Height;
                
                % Generate random x and y coordinates within the rectangle
                x_random = x_min + Width * rand(round(num_points), 1);
                y_random = y_min + Height * rand(round(num_points), 1);

                % Define grid edges based on your coordinate ranges
                x_edges = min(x_random):max(x_random)+1;
                y_edges = min(y_random):max(y_random)+1;

                data0_x = x_random+0.5;data0_y = y_random+0.5;
                                                                               
                % % Display the Random background image (debug only)
                % % Bin localizations into a grid
                % [N, ~, ~,~,~] = histcounts2(x_random, y_random, x_edges, y_edges);
                % figure; 
                % imshow(N,[],'InitialMagnification',1000);hold on;                
                % % scatter(data0_y,data0_x,5, 'filled')                
                % hROI0 = drawrectangle('Position',[Width/2-w/params.dxy/2,Height/2-h/params.dxy/2,w/params.dxy,h/params.dxy],...
                %     'LineWidth', 2, 'Color', 'r'); % Interactively draw an ROI
                % title('Random background');
                % 
                % % Get the mask from the drawn ROI
                % mask0 = createMask(hROI0); % Logical mask of the selected ROI  

                % Compute 2D histogram counts
                [N, ~, ~, ~, ~] = histcounts2(x_random, y_random, x_edges, y_edges);
                
                % Define the position and size of the rectangle directly
                roiPosition = [Width/2-w/params.dxy/2, Height/2-h/params.dxy/2, w/params.dxy, h/params.dxy];
                
                % Create a logical mask based on roiPosition
                mask0 = false(size(N)); % Initialize mask with zeros
                roiXRange = round(roiPosition(1)):(round(roiPosition(1)) + round(roiPosition(3)));
                roiYRange = round(roiPosition(2)):(round(roiPosition(2)) + round(roiPosition(4)));
                mask0(roiYRange, roiXRange) = true;

                % Define parameters
                params0.dr = 1;
                params0.calib = params.calib/params.dr*params0.dr;                
                params0.dxy = params0.dr;
                params0.num_pts = params.num_pts;
                              
                % Calculate autocorrelation using binned data inside the ROI
                [corrs0, params0] = ParCalculateXcorrApp(data0_x, data0_y, mask0, params0);
                
                
                % % Plot G(r) of this ROI (debug only)
                % figure('Units','inches','Position',[7.1562 3.8854 5.8333 2.7812]);
                % hold on;
                % errorbar(params.radii,corrs.C_11,corrs.dC_11,'s-','MarkerFaceColor','b','MarkerEdgeColor','b','Color','b','MarkerSize',2)
                % errorbar(params0.radii*params.dr,corrs0.C_11,corrs0.dC_11,'s--','MarkerFaceColor','k','MarkerEdgeColor','k','Color','k','MarkerSize',2)
                % plot([0 max(params.radii)], [1 1], 'k--', 'LineWidth', 1, 'HandleVisibility','off');
                % ylabel('G(r)');xlabel('r(nm)');
                % axis square
                
                if sum(isnan(corrs.C_11)) > 0 || sum(isnan(corrs0.C_11)) > 0
                    fprintf('Contain NaN, select ROI again. \n');
                else
                    g_r = [g_r; corrs.C_11];
                    g_r0 = [g_r0; corrs0.C_11];                              
                    mask_num = mask_num - 1;                    
                end
            end
        end
    end
    MeanGrScore = mean(g_r,1);
    MeanGr0Score = mean(g_r0,1);
    GrScore = g_r;
    Gr0Score = g_r0;
end

function ImOut_density = density_PALM(LocTable,LocTable_Xpos,LocTable_Ypos,loc_state,loc_XYLim)
    dxy = 30; 
    sigma_render = 30;
    Xmin = loc_XYLim(1);Xmax=loc_XYLim(2);Ymin = loc_XYLim(3);Ymax=loc_XYLim(4);
    Edges{1}=Xmin:dxy:Xmax+dxy; % make edge cover all data points
    Edges{2}=Ymin:dxy:Ymax+dxy;
    %% generate density based PALM image
    switch loc_state
        case 'all'
            Xpos = LocTable_Xpos;
            Ypos = LocTable_Ypos;
        case 'immobile'
            Xpos = LocTable_Xpos(LocTable.Area <= app.minMobileSRarea);
            Ypos = LocTable_Ypos(LocTable.Area <= app.minMobileSRarea);
    end
                
    Edges{1}=Xmin:dxy:Xmax+dxy; % make edge cover all data points
    Edges{2}=Ymin:dxy:Ymax+dxy;
                
    [Im,Xedges,Yedges,binX,binY] = histcounts2(Xpos,Ypos,Edges{1},Edges{2}); % bin every localizations into grid, similar to Im = hist3([Xpos',Ypos'],'Edges',Edges);

    TempX=-round(3*sigma_render/dxy)*dxy:dxy:round(3*sigma_render/dxy)*dxy;
    TempY=-round(3*sigma_render/dxy)*dxy:dxy:round(3*sigma_render/dxy)*dxy;
    
    ConVecX = exp(-0.5*(TempX/sigma_render).^2);
    ConVecX=ConVecX/sum(ConVecX);
    ConVecY = exp(-0.5*(TempY/sigma_render).^2);
    ConVecY=ConVecY/sum(ConVecY);
    % Cisse's verions of PALM rendering
%             Im2 = conv2(ConVecX,ConVecY,Im);
%             Im2=Im2/dxy/dxy*10^6; % convert unit of localizations in each bin from nm^-2 to um^-2
%             extra_pixels=(size(Im2)-size(Im))/2; % remove extra padding pixels
%             Im2=Im2((extra_pixels(1)+1):(end-extra_pixels(1)),(extra_pixels(2)+1):(end-extra_pixels(2)));
%             Im2=Im2(:,end:-1:1)';
%             if flipped
%                 ImOut_density = flip(Im2, 1);
%             else
%                 ImOut_density = Im2; 
%             end
    
    % ZW's updated shorten version, which gives same results
    Im2 = conv2(ConVecX,ConVecY,Im,'same');
    Im2=Im2/dxy/dxy*10^6; % convert unit of localizations in each bin from nm^-2 to um^-2
    ImOut_density = transpose(Im2); % x-y direction of histcounts2 are transposed with original image

end

function y = customModel(param, x)
% Define the model function
    % Ensure s is positive to avoid issues with negative exponentiation
    % s = abs(s);
    
    A = param(1); theta = param(2); ss= param(3); rho_p= param(4); 
    
    % Define the two exponential functions
    exp1 = A*exp(-x/theta)+1;
    exp2 = (1/(4*pi*ss))*exp(-x.^2 / (4*ss));
    
    
    % Convolve the two functions
    gr_psf = conv(exp1, exp2, 'same'); % Center result to match xData length

    gr_stoch = 1/(4*pi*ss*rho_p)*exp(-x.^2/(4*ss));

    y = gr_stoch+gr_psf;
    
    % % Scale the model data to fit the actual data range
    % yModel = convResult / max(convResult(:)); % Normalize if needed
end
