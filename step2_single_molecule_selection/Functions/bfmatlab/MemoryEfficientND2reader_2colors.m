function [ Color1Tiff,  Color2Tiff] = MemoryEfficientND2reader_2colors( nd2_file_to_open )
%MEMORYEFFICIENTND2READER 
%   MATLAB is annoying when it comes to parallel computing: cannot use the
%   BioFormat Memoizer wrapper in a parfor loop even though it works
%   perfectly well in a normal for loop. 
%   written by Anders Sejr Hansen

% add the neccesary paths:
addpath(genpath(['.' filesep 'Batch_MTT_code' filesep])); % MTT & BioFormats

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

% Define empty TIFF-like matrices. Very annoyingly, Nikon Elements make every other frame each color. So
% make two new tiff stacks for each color that are half the length in
% the 3rd dimension:
Color1Tiff = zeros(size(first_frame,1), size(first_frame,2), round(TotFrames/2));
Color2Tiff = zeros(size(first_frame,1), size(first_frame,2), round(TotFrames/2));

% now read in only one frame at a time to save on memory
for FrameIter = 1:size(Color1Tiff,3)
    Color1Tiff(:,:,FrameIter) = bfGetPlane(r,2*(FrameIter-1)+1);
    Color2Tiff(:,:,FrameIter) = bfGetPlane(r,2*(FrameIter-1)+2);
end

% FINISH
% Close the reader
r.close()

end

