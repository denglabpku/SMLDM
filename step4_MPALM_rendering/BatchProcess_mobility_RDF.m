% This script is to analyze MPALM mobility distribution plus membrane
% protein RDF.

% The input file is from BatchInput_mobility_RDF.m

%% Run mobility distribution fitting
close all;
D_SR_const = 1162 ; % The constant value from PT-area vs ground-truth D calibration curve

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
    
%% Auxiliary functions

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

