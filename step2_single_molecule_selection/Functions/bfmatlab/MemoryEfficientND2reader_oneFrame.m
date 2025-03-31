function [imgs_2d_matrix, TotFrames] = MemoryEfficientND2reader_oneFrame(nd2_file_to_open, varargin)
    % MemoryEfficientND2reader_oneFrame
    % Read desired frame from a large single-channel image stack. If no frame indx (varargin is empty)
    % is provided, then use the last frame.
    
    %MEMORYEFFICIENTND2READER 
    %   MATLAB is annoying when it comes to parallel computing: cannot use the
    %   BioFormat Memoizer wrapper in a parfor loop even though it works
    %   perfectly well in a normal for loop. 
    %   written by Anders Sejr Hansen
    
    
    % add the neccesary paths:
    % addpath(genpath(['.' filesep 'Batch_MTT_code' filesep])); % MTT & BioFormats
    
    %%% read ND2-file:
    
    %%%%%%%%%%%% MODIFY BY ZUHUI BUG FIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Disable chunkmap to read small ND2 correctly, temp solution until bug fix release
    autoloadBioFormats = 1;
    status = bfCheckJavaPath(autoloadBioFormats);
    assert(status, ['Missing Bio-Formats library. Either add bioformats_package.jar '...
        'to the static Java path or add it to the Matlab path.']);
    options = loci.formats.in.DynamicMetadataOptions();
    options.set("nativend2.chunkmap", "false");
    r = bfGetReader();
    r.setMetadataOptions(options);
    r = loci.formats.Memoizer(r);
    r.setId(nd2_file_to_open);
    TotFrames = r.getImageCount();
    %%%%%%%%%%%% MODIFY BY ZUHUI BUG FIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % % Construct an empty Bio-Formats reader
    % r = bfGetReader();
    % % Decorate the reader with the Memoizer wrapper
    % r = loci.formats.Memoizer(r);
    % % Initialize the reader with an input file
    % % If the call is longer than a minimal time, the initialized reader will
    % % be cached in a file under the same directory as the initial file
    % % name .large_file.bfmemo
    % r.setId(nd2_file_to_open);
    
    if isempty(varargin)
        % First figure out how many frames:
        % TotFrames = r.getImageCount();
        frame_idx = TotFrames;
    else 
        frame_idx = varargin{1};
    end
    
    % Read in the first frame to get size of images:
    imgs_2d_matrix = bfGetPlane(r, frame_idx);
    
    % FINISH
    % Close the reader
    r.close()
    
    end
    
    