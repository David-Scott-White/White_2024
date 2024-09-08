function [path, videos, time_s] = loadMMImages()
% For the case when the images are saved seperately, instead of as a image
% stack in MicroManger.
% Folder should be name "Default"
%  > all images as seperate .tif files
%  > metadata .txt file with the
%
% David S. White 
% 2021-09-20
% -------------------------------------------------------------------------

path = uigetdir('Select a Directory'); % make option to just give a path
listing = dir(path);

% Read metadata file
disp('>> Loading metadata.txt...');
fileID = fopen([listing(1).folder, '/', 'metadata.txt'], 'r');
metaData = textscan(fileID,'%s');
metaData = metaData{1};
fclose(fileID);

% Get number of channels
disp('>> Gathering info...')
numChannels = [];
numFrames = [];
width = [];
height = [];
for i = 1:length(metaData)
    if isempty(numChannels)
        if contains(metaData{i}, '"Channels":')
            numChannels = str2double(metaData{i+1}(1:end-1));
        end
    end
    if isempty(numFrames)
        if contains(metaData{i}, '"Frames":')
            numFrames = str2double(metaData{i+1}(1:end-1));
        end
    end
    if ~isempty(numChannels) && ~isempty(numFrames)
        break
    end
end
totalFrames = numFrames*numChannels;

% Get time
x = 1; 
elapsedTime = zeros(totalFrames,1); 
for i = 1:length(metaData)
    if contains(metaData{i}, '"ElapsedTime-ms":')
        elapsedTime(x,1) = str2double(metaData{i+1}(1:end-1));
        x = x +1; 
    end
end

% get videos
allVideos = [];
h = waitbar(0, 'Loading video...');
x = 1; 
for i = 1:length(listing)
    if contains(listing(i).name, '.tif')
        image = loadtiff([listing(i).folder,'/', listing(i).name]);
        if isempty(allVideos)
            allVideos = zeros(size(image,1), size(image,2), totalFrames); 
        end
        allVideos(:,:,x) = image + 1; % TEMP PATCH FOR BACKGROUND SUBTRACTION
        x = x+1; 
    end
    waitbar(i/totalFrames, h);
end
close(h)

% parse by channels
if numChannels > 1
    disp('>> Parsing videos by channel...')
    channelIdx = repelem(1:numChannels, 1, numFrames);
    videos = cell(numChannels,1);
    time_s = cell(numChannels,1);
    for i = 1:numChannels
        videos{i} = allVideos(:,:,channelIdx == i);
        time_s{i} = elapsedTime(channelIdx == i);
    end
end

% reformat time in [N,3 matrix] of [frame, real time (ms), elapsed time (s)]
for i = 1:numChannels
    temp = time_s{i}; 
    timeMat = zeros(length(temp), 3);
    timeMat(:,1) = 1:1:length(timeMat);
    timeMat(:,2) = temp;
    timeMat(:,3) = (temp-temp(1))/1000; % ms to seconds
    time_s{i} = timeMat; 
end
disp('>> Video loaded.')

% Channel Selection
videosTemp = videos; 
timeTemp = time_s;
videos = {};
time_s = {};
channelIdx = selectChannelWindow(numChannels);
if channelIdx(1) <= numChannels && channelIdx(2) <= numChannels
    videos{1} = videosTemp{channelIdx(1)};
    videos{2} = videosTemp{channelIdx(2)};
    time_s{1} = timeTemp{channelIdx(1)};
    time_s{2} = timeTemp{channelIdx(2)};
elseif channelIdx(1) > numChannels && channelIdx(2) <= numChannels
    videos{1} = videosTemp{channelIdx(2)};
    time_s{1} = timeTemp{channelIdx(2)};
elseif channelIdx(2) > numChannels && channelIdx(1) <= numChannels
    videos{1} = videosTemp{channelIdx(1)};
    time_s{1} = timeTemp{channelIdx(1)};
else
    msgbox('No Channels Selected');
end
end


