function [images, time_s] = glimpseToTif(path, boundingBox)
% -------------------------------------------------------------------------
% Author: David S. White
% Contact: dswhite2012@gmail.com
% License: GPLv3.0
% 
% Convert .glimpse video to .tif. 
%
% Updates: 
%   2021-06-04  DSW wrote the code
%   2021-06-07  DSW update to extract time (s)
% 
% Arguments: 
%   path = file path of video to load
%   splitImage = boolean. Crops images **hard coded from CoSMoS 1.0
% 
% Output: 
%   images = cell of K by 1 channels. loadedVideos{K} = [x,y ,N]
%   time_s = [3,N] (frame, time_physical_ms, time_s)
% 
% Note: 
%   Also saves a .tif for "all_images.tif" and seperate 
%   "all_images_green.tif", "all_images_red.tif" if splitImage == 1
% -------------------------------------------------------------------------

% --- Check Arguments
if nargin<1 || isempty(path)
    path = uigetdir();
end
if nargin < 2 || isempty(boundingBox)
    boundingBox = [];
end

% saveChannel = 'both' % not yet function (green, red, both)

% --- Load Glimpse File
foldstruc.gfolder = [path, '/'];
loadedGlimpseFile = load([path,'/header.mat']);
vid = loadedGlimpseFile.vid; 
disp(['>> ', path]);

% --- Possible size error, check and display warning. (10K max standard)
% ------ will need to write a correction function to save... 
if vid.nframes > 6500
    disp('Warning: Max Frames exceeds. Cut off at 6500 [fix in later version]'); 
    nFrames = 6500; 
else
    nFrames = vid.nframes;
end

% --- Process the glimpse file
images = zeros(vid.width,vid.height, nFrames);
time_s = zeros(nFrames, 3); 
for i = 1:nFrames
    fid=fopen([foldstruc.gfolder, num2str(vid.filenumber(i)),'.glimpse'],'r','b');
    fseek(fid,vid.offset(i),'bof');
    pc=fread(fid,[vid.width,vid.height],'int16=>int16');
    images(:,:,i) = uint16(pc+32768);  % <- int16 to unint16
    fclose(fid);
end
fclose('all');
time_s(:,1) = 1:nFrames; 
time_s(:,2) = vid.ttb(1:nFrames)'; 
time_s(:,3) = (time_s(:,2) - time_s(1,2))/ 1e3; % ms to seconds

% --- Save Images and Time (will not overwrite file if it already exists...)
try
    saveastiff(uint16(images),[foldstruc.gfolder,'all_images.tif']);
catch ME
    switch ME.message
        case 'File already exists.'
            disp('>> Warning: File already exists')
    end
end
save([path,'/time_s.mat'], 'time_s');

% ---  Auto crop Dual View EMCCD
if ~isempty(boundingBox)
    images_all = images;
    images = cell(2,1); 
    
    disp('>> Splitting Images...')
    
    fileNames = cell(2,1);
    fileNames{1} = 'all_images_green.tif';
    fileNames{2} = 'all_images_red.tif';
    
    % hard coded for CoSMoS 1.0. Values relative to array
    % [y start (from top to bottom), x start (left to right), imHeight, imWidth]
    
    %boundbox1 = [25, 70, 215, 215];
    %boundbox2 = [55 295, 215, 215];
    %boundingBoxes = [boundbox1; boundbox2];
    
    [imHeight, imWidth, nFrames] = size(images_all);
    if mod(imHeight,2) 
        imHeight = imHeight-1; 
    end
    if mod(imWidth, 2) 
        imWidth = imWidth -1;   
    end
    
    % [y start (from top to bottom), x start (left to right), imHeight, imWidth]
    bb = zeros(2,4);
    if imWidth > imHeight
        % images are horizontal
        % cropField = [row start, row stop, column start, column stop]
        bb(1,:) = [1, imHeight, 1, imWidth/2];
        bb(2,:) = [1, imHeight, imWidth/2+1, imWidth];
    end
    images = cell(2,1); 
    for i = 1:2
        images{i} = zeros(imHeight, imWidth/2, nFrames);
    end
    % vertical images..
    
    for i = 1:2
       
        r0 = bb(i,1); 
        r1 = bb(i,2); 
        c0 = bb(i,3); 
        c1 = bb(i,4);
        
        % Crop and Save
        for j = 1:nFrames
            images{i}(:, :, j) = images_all(r0:r1, c0:c1, j);
        end
        try
            saveastiff(uint16(images{i}),[foldstruc.gfolder,fileNames{i}]);
        catch ME
            switch ME.message
                case 'File already exists.'
                    disp('>> Warning: File already exists')
            end
        end
    end
end
end