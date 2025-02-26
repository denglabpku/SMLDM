function z = gaussian2D(par,xy)
% xy(:,:,1) = xData;
% xy(:,:,2) = yData;

% compute 2D gaussian
oneM = ones(1,size(xy,2));
z = par(:,7) + ...
    par(:,1).*exp(-(((xy(:,:,1)-par(:,5)*oneM).*cosd(par(:,2)*oneM)+(xy(:,:,2)-par(:,6)*oneM).*sind(par(:,2)*oneM))./(par(:,3)*oneM)).^2-...
    ((-(xy(:,:,1)-par(:,5)*oneM).*sind(par(:,2)*oneM)+(xy(:,:,2)-par(:,6)*oneM).*cosd(par(:,2)*oneM))./(par(:,4)*oneM)).^2);

% % compute 2D gaussian
% z = par(7) + ...
%     par(1)*exp(-(((xy(:,:,1)-par(:,5)).*cosd(par(2))+(xy(:,:,2)-par(6)).*sind(par(2)))./par(3)).^2-...
%     ((-(xy(:,:,1)-par(5)).*sind(par(2))+(xy(:,:,2)-par(6)).*cosd(par(2)))./par(4)).^2);

end

