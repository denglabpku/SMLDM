close all;clear;clc;
root_path = 'H:\2024\20240623_U2OS_FOXA2Halo-H2BmEosEM_dualColor_MPALM';
img_path = '20240623_Clust01_FOXA2-H2B_25nMPA646_30p5ms_10kframe_01.nd2';

image_file = fullfile(root_path,img_path);
frame_id = 1;
[Ch1_img,Ch2_img] = bfopenframe_twoChannel(image_file,frame_id);

%% 560 channel PALM

x_Ch1_beads = [];
y_Ch1_beads = [];

figure('Name','First frame of Channel 1');
ax = axes('Position',[0 0 1 1]);
imshow(Ch1_img,[],'Parent',ax);
h1 = imcontrast();
uiwait(h1);
seg_prompt = {'Enter number of beads:'};
seg_prompt_name = 'Wait for user';
numlines = 1;
seg_defaults = {'1'};  % default is one nucleus mask, and zero compartments mask
options.Resize = 'on';
options.WindowStyle = 'normal';
seg_answer = inputdlg(seg_prompt,seg_prompt_name,numlines,seg_defaults);
beads_number = str2double(seg_answer{1});  % number of nucleus mask


for i = 1:beads_number
    roi = drawrectangle(ax,"FaceAlpha",0,"Label",num2str(i));
    round_pos = round(roi.Position);
    crop_img = Ch1_img(round_pos(2):round_pos(2)+round_pos(4),round_pos(1):round_pos(1)+round_pos(3));
    
    [xo,yo,zo] = prepareSurfaceData(round_pos(1):round_pos(1)+round_pos(3),round_pos(2):round_pos(2)+round_pos(4),crop_img); 

    fitresult = fit_gaussian2(xo,yo,zo);

    x_Ch1_beads(i) = fitresult(4);
    y_Ch1_beads(i) = fitresult(5);


end
hold(ax,"on")
scatter(ax,x_Ch1_beads,y_Ch1_beads,'o');hold(ax,'off')

xy_Ch1_beads = [x_Ch1_beads',y_Ch1_beads'];

%% 642 channel MPALM

x_Ch2_beads = [];
y_Ch2_beads = [];

figure('Name','First frame of Channel 2');
ax = axes('Position',[0 0 1 1]);
imshow(Ch2_img,[],'Parent',ax);
h1 = imcontrast();
uiwait(h1);

for i = 1:beads_number
    roi = drawrectangle(ax,"FaceAlpha",0,"Label",num2str(i));
    round_pos = round(roi.Position);
    crop_img = Ch2_img(round_pos(2):round_pos(2)+round_pos(4),round_pos(1):round_pos(1)+round_pos(3));
    
    [xo,yo,zo] = prepareSurfaceData(round_pos(1):round_pos(1)+round_pos(3),round_pos(2):round_pos(2)+round_pos(4),crop_img); 

    fitresult = fit_gaussian2(xo,yo,zo);

    x_Ch2_beads(i) = fitresult(4);
    y_Ch2_beads(i) = fitresult(5);

end
hold(ax,"on")
scatter(ax,x_Ch2_beads,y_Ch2_beads,'o');hold(ax,'off')

xy_Ch2_beads = [x_Ch2_beads',y_Ch2_beads'];


%% Align two channel
MOVING = Ch1_img;
FIXED = Ch2_img;

mp = xy_Ch1_beads;
fp = xy_Ch2_beads;

t = fitgeotform2d(mp,fp,"similarity");    

Rfixed = imref2d(size(FIXED));
registered = imwarp(MOVING,t,OutputView=Rfixed);

% figure
% imshowpair(FIXED,registered,"blend")

progressbarText(0);
for i = 1:total_frame
    

    tform = simtform2d(Scale(i),RotationAngle(i),[Translation_1(i) Translation_2(i)]);
    row_idx = app.mainapp.Ch2_loc_table.Frame == i;

    [u,v] = ...
        transformPointsForward(tform,app.mainapp.Ch2_loc_table.xPix(row_idx),app.mainapp.Ch2_loc_table.yPix(row_idx));

    app.mainapp.Ch2_loc_table.regis_Xpos(row_idx) = u*app.mainapp.camPixelSize_Ch1.Value;                
    app.mainapp.Ch2_loc_table.regis_Ypos(row_idx) = v*app.mainapp.camPixelSize_Ch1.Value;
    progressbarText(i/total_frame);                               
end
disp('test!')






%%
function fitresult = fit_gaussian2(xData,yData,zData)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% assign start point
EmissionWavelength = 664;
impars.PixelSize=0.11; % um per pixel, Andor iXon Ultra 897
impars.psf_scale=1.35; % PSF scaling
impars.wvlnth= EmissionWavelength/1000; %emission wavelength in um
impars.NA=1.49; % NA of detection objective
impars.psfStd= impars.psf_scale*0.55*(impars.wvlnth)/impars.NA/1.17/impars.PixelSize/2; % PSF standard deviation in pixels

% assign the fit model using MTT paper
% g = fittype('(I/r0) * exp(-(1/(2*r0^2))*((y - y0)^2 + (x - x0)^2)) + 498','independent',{'x','y'},'dependent','z');
g = fittype('I * exp(-(1/(2*r0^2))*((y - y0)^2 + (x - x0)^2)) + m','independent',{'x','y'},'dependent','z');
fo = fitoptions(g);
fo.StartPoint = [max(zData)-min(zData),min(zData),impars.psfStd,mean(xData),mean(yData)];
% fo.StartPoint = [max(zData)-498,impars.psfStd,mean(xData),mean(yData)];
fo.Lower = [0,0,0,min(xData),min(yData)];
fo.Upper = [max(zData),mean(zData),Inf,max(xData),max(yData)];
% fo.Lower = [min(zData),0,0,0];
% fo.Upper = [max(zData),impars.psfStd+1,512,512];
fitobject  = fit([xData,yData],zData,g,fo);
fitresult = [fitobject.I, fitobject.r0, fitobject.m,fitobject.x0,fitobject.y0];
% fitresult = [fitobject.I, fitobject.r0,fitobject.x0,fitobject.y0];

% % % Check the fit result
% figure
% plot(fitobject, [xData yData], zData);
% xlabel('X-coord (pixel)');ylabel('Y-coord (pixel)');zlabel('Intensity (A.U.)');
end


