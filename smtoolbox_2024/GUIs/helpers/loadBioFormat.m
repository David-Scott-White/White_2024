function [images, time_s] = loadBioFormat(file, path, splitChannels, saveImages)
% -------------------------------------------------------------------------
% Author: David S. White
% Contact: dswhite2012@gmail.com
% License: GPLv3.0
%
% Wrapper around bioformats to load images from micro-manager (.ome.tif).
%   see: https://docs.openmicroscopy.org/bio-formats/6.1.0/users/matlab/index.html
%
% Updates:
%   2021-06-07  DSW wrote the code
%
% Requirements:
%   bfmatlab
%
% Arguments:
%   filePath = file path of video to load
%   splitChannels = boolean. [not yet functional]
%
% Output:
%   images = cell of K by 1 channels. loadedVideos{K} = [x,y ,N]
%   time_s = cell of K by 1 channels. time_s{K} = [3,N] (frame, time_physical_ms, time_s)
% -------------------------------------------------------------------------

% Check Arguments
if nargin < 1 || isempty(file)
    [file, path] = uigetfile('*.tif');
    filePath = [path,file];
    splitChannels = 1;
    saveImages = 0;
    
end
if nargin < 3
    splitChannels = 1;
    saveImages = 1;
end
filePath = [path,file];

% option for merge two?

% Load image (turn off warnings?)
reader = bfGetReader(filePath);
omeMeta = reader.getMetadataStore();
data = bfopen(filePath);

% To do- add a waitbar for processing images after loaded.

% get pixel size
imageData = data{1};
numFrames = size(imageData, 1);
[x,y] = size(imageData{1});

% Check for multiple channels          ***(hard code max of 2 for now...)
channelIdx = zeros(numFrames, 2);
channelIdx(:,1) = 1:numFrames;
channelCheck = strfind(imageData{1,2}, 'C=');
if isempty(channelCheck)
    numChannels = 1;
    channelIdx(:,2) = 1;
else
    % reported as 'C=1/X) -> edit into a better loop using numChannels
    numChannels = str2double(imageData{1,2}(channelCheck+4));
    if numChannels == 2
        channelIdx(1:2:numFrames-1, 2) = 1;
        channelIdx(2:2:numFrames, 2) = 2;
    elseif numChannels == 4
        channelIdx(1:4:numFrames-1, 2) = 1;
        channelIdx(2:4:numFrames, 2) = 2;
        channelIdx(3:4:numFrames,2) = 3;
        channelIdx(4:4:numFrames,2) = 4;
    end
end

% loop through all images.
stopIdx = 0;
images = cell(numChannels,1);
time_s = cell(numChannels,1);  % frame, real time (s), time_s - time_s 0
for k = 1:numChannels
    idx = channelIdx(channelIdx(:,2) == k, 1);
    images{k} = zeros(x,y,length(idx));
    time_s{k} = zeros(length(idx), 3);
    for i = 1:length(idx)
        images{k}(:,:,i) = imageData{idx(i)};
        time_s{k}(i,1) = i;
        
        % possible to lose metadata if the recording was manually stopped.
        try
            time_s{k}(i,2) = double(omeMeta.getPlaneDeltaT(0,idx(i)-1).value());
        catch ME
            switch ME.identifier
                case 'MATLAB:Java:GenericException'
                    disp(['>> No Metadata for image ', num2str(idx(i)), '. Removing...']);
                    stopIdx = i;
            end
        end
        if stopIdx < 0 % <- there is a bug here for an infinite loops somehow...
            images{k} =  images{k}(:,:,1:stopIdx-1);
            time_s{k} = time_s{k}(1:stopIdx-1, :);
        end
    end
end

% Grab first time and normalize all other times. Correct other times.
time_s_0 = time_s{1}(1,2);
for k = 1:numChannels
    time_s{k}(:,3) = (time_s{k}(:,2) - time_s_0)/1e3; % ms to seconds
    
    % clean up, might need a correction/optinos if distributing 
    dx = mean(diff(time_s{k}(:,3)));
    time_s{k}(:,3) = time_s{k}(:,3) + dx; % correct first time (non-zero)
end
disp(['loadBioFormat:  ', fileparts(filePath), ' loaded.']);

% temp, write as seperate function?
if saveImages
    for i = 1:numChannels
        saveastiff(uint16(images{i}),[path, file(1:end-8),'_C',num2str(i),'.tif']);
    end
end
