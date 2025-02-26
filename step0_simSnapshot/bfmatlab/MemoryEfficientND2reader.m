function [ imgs_3d_matrix ] = MemoryEfficientND2reader( nd2_file_to_open )
%MEMORYEFFICIENTND2READER 
%   MATLAB is annoying when it comes to parallel computing: cannot use the
%   BioFormat Memoizer wrapper in a parfor loop even though it works
%   perfectly well in a normal for loop. 
%   written by Anders Sejr Hansen

% add the neccesary paths:
% addpath(genpath(['.' filesep 'Batch_MTT_code' filesep])); % MTT & BioFormats

%%% read ND2-file:

% Construct an empty Bio-Formats reader
r = bfGetReader();
% Decorate the reader with the Memoizer wrapper
r = loci.formats.Memoizer(r);
% Initialize the reader with an input file
% If the call is longer than a minimal time, the initialized reader will
% be cached in a file under the same directory as the initial file
% name .large_file.bfmemo
r.setId(nd2_file_to_open);

% First figure out how many frames:
TotFrames = r.getImageCount();

% Read in the first frame to get size of images:
first_frame = bfGetPlane(r, 1);

% pre-allocate memory for the large 3D matrix
imgs_3d_matrix = zeros(size(first_frame,1), size(first_frame,2), round(TotFrames));

% go through all of the frames and read them in one at a time
for FrameIter = 1:TotFrames
    imgs_3d_matrix(:,:,FrameIter) = bfGetPlane(r,FrameIter);
end

% FINISH
% Close the reader
r.close()

end

