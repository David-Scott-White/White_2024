function rois = projectAOI(rois, video)
% Replaces projectTimeSeries 
% Use if there is no drift correction to the centroids, only video or none
%
% Input:
%   rois. stuct. see findAOI.m 
%   video. cell of n,1 where video{1} = [x,y,frames];
%
% Author: David S. White 
% Date 2022-02-25
% License: MIT
% -------------------------------------------------------------------------

[nROIs, nChannels] = size(rois);
nTotal = nChannels*nROIs;
p = 1;
wb = waitbar(0, 'Projecting AOIs...');
removeIdx = [];
for j = 1:nChannels
    videoIter = video{j};
    nFrames = size(videoIter, 3);
    for i = 1:nROIs
        % grab time series (integrate over the AOI)
        rois(i,j).timeSeries = zeros(nFrames,1);
        pixelList = rois(i,j).pixelList;
        numPixels = size(pixelList,1);
        for k = 1:numPixels
            col = rois(i,j).pixelList(k,1);
            row = rois(i,j).pixelList(k,2);
            rois(i,j).timeSeries = rois(i,j).timeSeries + double(squeeze(video{j}(row,col,:)));
            
            % zero
             rois(i,j).timeSeries = rois(i,j).timeSeries-min(rois(i,j).timeSeries);
        end

        % store the spot (add modSize to AOI area)
        modArea = 2; % if [3,3], 1 makes 5x5. 2 makes 7x7 
        bb = rois(i,j).boundingBox(1,:); 
        bb(1) = bb(1)-modArea; 
        bb(2) = bb(2)-modArea;
        bb(3:4) = bb(3:4)+modArea*2; 
        pixelList = boundBoxToPixels(bb); % need var option to control for area size
        nPixelList = size(pixelList,1)^0.5; 
        pixelValues = zeros(length(pixelList), nFrames);
        if sum(sum(pixelList<=0)) || sum(sum(pixelList> 512))
            removeIdx = [removeIdx; i];
        else
            for c = 1:length(pixelList)
                pixelValues(c,:) = videoIter(pixelList(c,2), pixelList(c,1), :);
            end
            rois(i,j).spot = reshape(pixelValues,[nPixelList, nPixelList, nFrames]);
        end
        waitbar(p/nTotal, wb);
        p = p +1;
    end
end
rois(removeIdx, :) = []; 
close(wb)

[nROIs, nChannels] = size(rois)

% need some way to detect beads

% Remove ROIs that drop to 0 (drift correction)
% rois(removeIdx,:) = [];