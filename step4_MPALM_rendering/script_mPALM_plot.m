% testA = smoothedA;
% testA(testA < 40 & testA > 0 ) = 40;
% 
% hF = figure();
% hI = imagesc(testA);
% axis image;
% axis off;
% % Create custom colormap
% J = customcolormap([0 0.5 1], [0 1 0; 1 0 0; 0 0 1]);
% J(1,:) = [0 0 0];
% colormap(J);
% % Define values for colormap mapping
% values = [39 60 80]; % Blue, green, red thresholds
% % Render image using imagesc and set color axis limits
% clim([min(values), max(values)]);
% ax = gca;
% ax.Position = [0 0 1 1]
% 
% % Convert figure into RGB image
% frame = getframe(hF);
% rgbImage = frame2im(frame);
% 
% C = imfuse(nucleus_partition,rgbImage,'blend','Scaling','joint');
% figure;imshow(C)
% 
% figure;imagesc(nucleus_partition); axis image

%%
testA = smoothedA;
testA(testA < 40 & testA > 0 ) = 40;
hF = figure();
set(hF, 'Units', 'pixels');
figSize = size(testA);
figWidth = figSize(2);
figHeight = figSize(1);
set(hF, 'Position', [0 0 figWidth figHeight]);

hI = imagesc(testA);
axis image;
axis off;

% Create custom colormap
J = customcolormap([0 0.5 1], [0 1 0; 1 0 0; 0 0 1]);
J(1,:) = [0 0 0];
colormap(J);

% Define values for colormap mapping
values = [39 60 80]; % Blue, green, red thresholds

% Render image using imagesc and set color axis limits
caxis([min(values), max(values)]);

ax = gca;
ax.Position = [0 0 1 1]

% Convert figure into RGB image
frame = getframe(hF);
rgbImage = frame2im(frame);

% Resize RGB image to the same size as testA
rgbImage = imresize(rgbImage, [figHeight figWidth]);

figure
% Display the RGB image using imshow
imshow(nucleus_partition);

figure;imshow(rgbImage)

C = imfuse(nucleus_partition,rgbImage,'blend','Scaling','joint');
figure;imshow(C)

figure;imshow(nucleus_partition);freezeColors;
mobile_raw_table = raw_table(raw_table.Area > 50,:);
Xpos_rescale = (mobile_raw_table.Xpos-RenderParam.Xrange(1))./RenderParam.dxy+0.5 ;
Ypos_rescale = (mobile_raw_table.Ypos-RenderParam.Yrange(1))./RenderParam.dxy+0.5 ;
Area = mobile_raw_table.Area;
% overlay detection localization
hold on;
scatter(Ypos_rescale,Xpos_rescale,15,Area,'filled');
J = customcolormap([0 0.5 1], [0 1 0; 1 0 0; 0 0 1]);
colormap(J);
% Define values for colormap mapping
values = [40 60 80]; % Blue, green, red thresholds
% Render image using imagesc and set color axis limits
clim([min(values), max(values)]);


figure;imshow(C);freezeColors;
mobile_raw_table = raw_table(raw_table.Area > 0,:);
Xpos_rescale = (mobile_raw_table.Xpos-RenderParam.Xrange(1))./RenderParam.dxy+0.5 ;
Ypos_rescale = (mobile_raw_table.Ypos-RenderParam.Yrange(1))./RenderParam.dxy+0.5 ;
Area = mobile_raw_table.Area;
% overlay detection localization
hold on;
scatter(Ypos_rescale,Xpos_rescale,15,Area,'filled');
J = customcolormap([0 0.5 1], [0 1 0; 1 0 0; 0 0 1]);
colormap(J);
% Define values for colormap mapping
values = [40 60 80]; % Blue, green, red thresholds
% Render image using imagesc and set color axis limits
clim([min(values), max(values)]);

figure;
montage({nucleus_partition,rgbImage})
%%




















% Convert figure into RGB image
frame = getframe(hF);
rgbImage = frame2im(frame);

figure
% Display the RGB image using imshow
imshow(rgbImage);
