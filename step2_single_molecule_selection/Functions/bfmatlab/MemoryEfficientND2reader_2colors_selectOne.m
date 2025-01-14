function ColorTiff = MemoryEfficientND2reader_2colors_selectOne( nd2_file_to_open ,channel_order)
%MEMORYEFFICIENTND2READER 
%   MATLAB is annoying when it comes to parallel computing: cannot use the
%   BioFormat Memoizer wrapper in a parfor loop even though it works
%   perfectly well in a normal for loop. 
%   written by Anders Sejr Hansen
%     ZW modified to only output one channel of dual channel ND2 images
%     channel_order: best accurate way to select channel_order is to use
%     NIS-Viewer to check which channel appear first, and that channel will be
%     order 1, and so on.

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
ColorTiff = zeros(size(first_frame,1), size(first_frame,2), round(TotFrames/2));

% now read in only one frame at a time to save on memory
for FrameIter = 1:size(ColorTiff,3)
    ColorTiff(:,:,FrameIter) = bfGetPlane(r,2*(FrameIter-1)+channel_order);
end

% FINISH
% Close the reader
r.close()

end

