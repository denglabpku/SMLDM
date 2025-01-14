function [m_data1] = apply_mask(data1,mask,k)

if isempty(data1) %returns empty if no data 
    m_data1 = [];
return
end

X1 = data1(:,1); % x positions of data
Y1 = data1(:,2); % y positions of data

inders = 1:numel(X1); %numbered indexes of original data

NN = X1<size(mask,1); % getting rid of out of range indexes
WW = X1>0;
NM = Y1<size(mask,2); 
WA = Y1>0;

RERE = NN.*WW.*NM.*WA;
loggg1 = logical(RERE); %logical indexes of data range of mask
        
if sum(loggg1) == 0 %if there is no data, return everything empty
    m_data1 = [];     
return
end        
            
subs = [ceil(double(X1(loggg1))) ceil(double(Y1(loggg1)))];
% the x and y positions that are within the range of the mask
% subs has removed particles that are out of range of the mask      
% subs is ceiling in order to index directly into the mask 

struct_inds = inders(loggg1);
%the numerical indexes of the original data in range of mask

ch1_inds = 1 + (subs(:,1)-1)*k(1) + (subs(:,2)-1)*k(2);
% getting 1D index from 2D subscripts subs [x y] using a 
% workaround for sub2ind. sub2ind calls led to decreased performance        

keep1 = mask(ch1_inds); % applying the logical mask to the data     
inds_zeros = struct_inds'.*keep1; % indexes of all data within mask
% with zeros at indexes of data outside of mask        

final_inds = inds_zeros(inds_zeros>0);        
m_data1 = data1(final_inds,:); %the masked data1                      
end

