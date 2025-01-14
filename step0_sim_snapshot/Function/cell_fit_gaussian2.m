function fitresult = cell_fit_gaussian2(xData,yData,zData)
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
g = fittype('I * exp(-(1/(2*r0^2))*((y - y0)^2 + (x - x0)^2)) + 498','independent',{'x','y'},'dependent','z');
fo = fitoptions(g);
% fo.StartPoint = [max(zData)-min(zData),min(zData),impars.psfStd,mean(xData),mean(yData)];
fo.StartPoint = [max(zData)-498,impars.psfStd,mean(xData),mean(yData)];
% fo.Lower = [min(zData),0,0,0,0];
% fo.Upper = [max(zData),max(zData),impars.psfStd+1,512,512];
fo.Lower = [min(zData),0,0,0];
fo.Upper = [max(zData),impars.psfStd+1,512,512];
fitobject  = fit([xData,yData],zData,g,fo);
% fitresult = [fitobject.I, fitobject.r0, fitobject.m,fitobject.x0,fitobject.y0];
fitresult = [fitobject.I, fitobject.r0,fitobject.x0,fitobject.y0];

% % Check the fit result
% figure
% plot(fitobject, [xData1 yData1], zData1);
% xlabel('X-coord (pixel)');ylabel('Y-coord (pixel)');zlabel('Intensity (A.U.)');
end

