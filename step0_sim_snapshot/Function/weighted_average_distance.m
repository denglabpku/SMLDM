function [weighted_average_distance] = weighted_average_distance(localizations)
%weighted_average_distance Summary of this function goes here
%   Detailed explanation goes here

% TO-DO: modify this function to use photon numbers as weight

%========= version1: use pixel value as weight ========================== %
% Calculate the number of localizations
n = size(localizations, 1);

% Initialize a variable to store the sum of weighted distances
weighted_distance_sum = 0;
weight_sum = 0;

% Loop over all pairs of localizations
for i = 1:n
    for j = i+1:n
        % Calculate the distance between the i-th and j-th localizations
        distance = norm(localizations(i,1:2)-localizations(j,1:2));
        
        % Add the weighted distance to the sum
        weighted_distance_sum = weighted_distance_sum + localizations(i,3)*localizations(j,3)*distance;
        weight_sum = weight_sum + localizations(i,3)*localizations(j,3);
    end
end

% Finally, calculate the weighted average distance and display the result
weighted_average_distance = weighted_distance_sum / weight_sum;
% disp(weighted_average_distance)

%========= version2: use photon number as weight ========================== %
end

