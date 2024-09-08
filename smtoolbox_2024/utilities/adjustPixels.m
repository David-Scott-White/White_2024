function [pixels,boundingBox] = adjustPixels(centroid,radius)
% David S. white
% dswhite2012@gmail.com
% 2019-12-02

% compute pixels of bounding box around the center at the defined radius.
tlx = round(centroid(1)-radius); % top left x (was floor)
tly = round(centroid(2)-radius); % top left y % bounding box in region props format
px = tlx:(tlx+radius*2)-1; % pixles in x (subtract 1 for imcrop correction)
py = tly:(tly+radius*2)-1; % pixels in y 
pixels = []; 
for x =1:length(px)
    for y = 1:length(py)
        pixels = [pixels; [px(x),py(y)]];
    end
end
boundingBox = [tlx,tly,length(px)-1,length(py)-1];    

