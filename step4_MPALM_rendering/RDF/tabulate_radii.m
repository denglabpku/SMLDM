
function [histobam] = tabulate_radii(DATA1,DATA2,rad_dist_bins,histobam)

m_data1 = DATA1; %the masked data1
m_data2 = DATA2; %the masked data2  



for k = 1: numel(m_data1(:,1)) %cycling through only data points within mask

A = m_data1(k,:); %particle k's x and y position

B = A(ones(size(m_data2,1),1),:); %repeat particle 1's position 
%to match size of data2
diff = B - m_data2; %comparing particle k to all particles 
%in data2
diff_squared = diff.*diff; % [x_diff^2 y_diff^2]
r_diffs = (sqrt(diff_squared(:,1) + diff_squared(:,2)))'; 
%radial difference values


[new_hists, indexers2] = histc(r_diffs,rad_dist_bins); %binning by radius 

% rad_dist_bins = radial distribution 

%indexing the column indexes by k, (m_data1(k))
%amazingly, this index takes car of both channels indexing in an elegant
%way: all_radii_indexes(k,j) will return a vector containing which 
%bins the radii between m_data(k) and m_data(j) ((NOTE m_data is actually 
% N by 2 since it contains x and y values for position...)

histobam = histobam + new_hists; %adding the histograms for g(r)
  
end

end
