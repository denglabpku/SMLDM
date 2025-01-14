function plot_ksdensity_histogram_curve(data,color)
%ksdensity_histogram_curve Calculate smoothed histogram using ksdensity
%kernel
%   Detailed explanation goes here

% Compute the histogram
% [n,edges] = histcounts(data);

% % % Compute the kernel density estimate
% xgrid = linspace(min(data),max(data),100);
% % % xgrid = linspace(0,100,100);
% pd = fitdist(data,'Kernel','Width',0.1);
% y = pdf(pd,xgrid);

% Plot the histogram and the smoothed curve
% bar(edges(1:end-1),n,'hist')
% hold on
% plot(xgrid,y*numel(data)*(edges(2)-edges(1)),'Color',color,'LineWidth',2)
% plot(xgrid,y,'Color',color,'LineWidth',2)


% Estimate the PDF using kernel smoothing
[f, xi] = ksdensity(data,'Bandwidth',0.5);

% % Calculate the probability by multiplying the density values by bin widths
% binWidth = xi(2) - xi(1);
% probability = f * binWidth;
% fprintf('%g\n',sum(probability));

% Plot the kernel smoothed probability
% plot(xi, probability,'Color',color,'LineWidth',2);
plot(xi,f,'Color',color,'LineWidth',2)

% hold off
end

