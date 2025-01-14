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
[xvals, yvals] = meshgrid([-rmax:rmax],[-rmax:rmax]);

% transform to polar coordinates with v as image values 
[~,r,v] = cart2pol(xvals,yvals, zvals); % ~ = theta

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
for j = 1:rmax+1 %length(r),
    m = bin==j;
    n2 = sum(m);
    if n2==0, vals(j)=0; er(j)=0; 
    else
        vals(j) = sum(m.*vv)/n2;        
    end
end

end

