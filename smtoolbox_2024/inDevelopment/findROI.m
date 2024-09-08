function rois = findROI(image, info, showResult)
% Find regions of interset (ROI) (e.g., single molecules)
% Note, replaces findAOI for now. Just finds the AOI, does not filter
%
%
%
% Author: David S. White 
% Date: 2022-01-02
% License: MIT
% -------------------------------------------------------------------------

if nargin < 2 || isempty(info)
    info = struct;
    info.area = 4;
    info.radius = 2;
    info.fp = 15;
    info.gauss = 1;
    info.tol = 1e-6;
    info.method = 'GLRT';
    showResult = 0;
end

if nargin < 3
    showResult = 0;
end

switch info.method
    
    case 'GLRT'
        [BW, ~, FAR] = GLRTfiltering (image, info.area, info.radius, info.fp);
        BW = bwareaopen(BW, info.radius, 8);
        BW = bwareafilt(BW, [1, 10]);
        stats = regionprops(BW,FAR, 'WeightedCentroid','Area');
        
    case 'Otsu'
        se = strel('disk',3);
        im = imtophat(mat2gray(image), se);
        BW = imbinarize(im,'global');
        BW = bwareaopen(BW, info.radius, 8);
        stats = regionprops(BW,im,'WeightedCentroid','Area');
end

% Filter the ROI (need to set these up as inputs 
% drop anything less than area == 4 
dropAreaIdx = find(vertcat(stats.Area)<4);
stats(dropAreaIdx) = []; 

centroids = cat(1,stats.WeightedCentroid);
numROIs = size(centroids,1);

% Drop anything too close in distance
% 
spotsRemoved =[];
for i = 1:numROIs
    dist = sqrt((centroids(i,1) - centroids(:,1)).^2 + (centroids(i,2) - centroids(:,2)).^2);
    removeIndex = find(dist < info.radius*2);
    if length(removeIndex) > 1
        spotsRemoved = [spotsRemoved; removeIndex];
    end
end
if length(spotsRemoved) > 1
    disp(['Number of Spots Removed: ', num2str(length(spotsRemoved))])
end
if ~isempty(spotsRemoved)
    centroids(spotsRemoved,:) = [];
    numROIs = size(centroids,1);
end

% reformat bounding box
[boundingBox, pixeList] = makeBoundingBox(centroids, info.area);

% Gausian fit the ROIs





