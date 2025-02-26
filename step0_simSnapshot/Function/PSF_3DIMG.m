function PSF_3DIMG(i_plane, PSF_idx, cell_PSF)
%PSF_3DIMG Draw 3D raw PSF.
%   Detailed explanation goes here. 

%% Zuhui Wang
%% 2020/11/15
%%

d = cell_PSF.xyI{i_plane}{PSF_idx};
X = d(:,1); Y = d(:,2); Z = d(:,3);
Xs = unique(X);
Ys = unique(Y);
Xi = arrayfun(@(x) find(Xs == x), X); %return index of Xs that equal to X
Yi = arrayfun(@(y) find(Ys == y), Y); %return index of Ys that equal to Y
Li = Yi + (Xi -1) * numel(Ys);
XYZ = nan(numel(Ys), numel(Xs));
XYZ(Li) = Z;

% Plot the surface of PSF
fig = figure("Visible","on");
ax = gca();
surf(ax, Xs,Ys,XYZ);
% shading interp;
view(-37.5,52.8);
set(ax, 'YDir', 'reverse');
colorbar;
xlabel(ax, 'X coordinate (pixel)');
ylabel(ax, 'Y coordinate (pixel)');
zlabel(ax, 'Intensity (a.u.)');
title(ax, ['Extract ' num2str(PSF_idx) '-th PSF profile from ' num2str(i_plane) '-th image']); 

end

