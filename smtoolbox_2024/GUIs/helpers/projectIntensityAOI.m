function rois = projectIntensityAOI(rois, video)
% Project intensity of AOI over frames. assumes no change in the location 
% of the bounding box 


[nROIs, nChannels] = size(rois);
nTotal = nChannels*nROIs;
wb = waitbar(0, 'Projecting AOIs...');
p = 1;
for j = 1:nChannels
    [~, ~, nFrames] = size(video{j});
    for i = 1:nROIs
        rois(i,j).timeSeries = zeros(nFrames,1);
        pixelList = rois(i,j).pixelList;
        numPixels = size(pixelList,1);
        for k = 1:numPixels
            col = rois(i,j).pixelList(k,1);
            row = rois(i,j).pixelList(k,2);
            rois(i,j).timeSeries = rois(i,j).timeSeries + double(squeeze(video{j}(row,col,:)));
        end
        p = p +1;
        waitbar(p/nTotal, wb);
        p = p +1;
    end
end
close(wb)


end