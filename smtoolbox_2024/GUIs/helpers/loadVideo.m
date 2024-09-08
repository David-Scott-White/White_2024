function [file, path, videos, time_s] = loadVideo()
% David S. White 
% 2021-09-14

% General function to call for loading single molecule video from different
% sources (e.g, .tif, .ome.tif, .glimpse) 

videos = {}; 
time_s = {}; 

% need to include a tag to label video channels

[file, path] = uigetfile('*');
if contains(file,'.ome')
    splitChannels = 1; 
    saveImages = 0; 
    [allVideos, allTime] = loadBioFormat(file, path, splitChannels, saveImages); 
    numChannels = size(allVideos,1); 
    if numChannels > 1
        channelIdx = selectChannelWindow(numChannels);
        if channelIdx(1) <= numChannels && channelIdx(2) <= numChannels 
            videos{1} = allVideos{channelIdx(1)}; 
            videos{2} = allVideos{channelIdx(2)};
            time_s{1} = allTime{channelIdx(1)}; 
            time_s{2} = allTime{channelIdx(2)}; 
        elseif channelIdx(1) > numChannels && channelIdx(2) <= numChannels
            videos{1} = allVideos{channelIdx(2)};
            time_s{1} = allTime{channelIdx(2)};
        elseif channelIdx(2) > numChannels && channelIdx(1) <= numChannels
            videos{1} = allVideos{channelIdx(1)};
            time_s{1} = allTime{channelIdx(1)};
        else
            msgbox('No Channels Selected');
        end
    else
        videos = allVideos;
        time_s =allTime;
    end
else
    videos{1} = loadtiff([path, file]); 
    time_s{1} = 1:size(videos,3);
end

% Convert  to uint16 (or greyscale?)
for i = 1:length(videos)
    videos{i} = uint16(videos{i});
    
    if min(min(videos{i}(:,:,1))) == 0
        for j = 1:size(videos{i},3)
            videos{i}(:,:,j) = videos{i}(:,:,j) + 1;
        end
    end
end

end