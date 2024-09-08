function [rois, zeroIdx] = findAOI(image, aoi, showResult)
% David S. White
% 2021-09-19

% Function to find Areas of Interest (AOI) in an image.
% Code utilizes

if nargin < 2 || isempty(aoi)
    aoi = struct;
    aoi.area = 3;
    aoi.radius = 2;
    aoi.fp = 20;
    aoi.gauss = 1;
    aoi.tol = 1e-4;
    aoi.method = 'GLRT';
    aoi.minSeperation = aoi.radius;
end

if nargin < 3 || isempty(showResult)
    showResult = 0;
end

% MANUALLY SET ----------------------
% aoi.minSeperation = 5;
image = double(image);

switch aoi.method
   
    case 'GLRT'
        [BW, ~, FAR] = GLRTfiltering (image, aoi.area, aoi.radius, aoi.fp);
        BW = bwareaopen(BW, aoi.radius, 8);
        BW = bwareafilt(BW, [1, 10]);
        
        stats = regionprops(BW,FAR, 'WeightedCentroid', 'MeanIntensity');
        
    case 'Otsu'
        se = strel('disk',3);
        im = imtophat(mat2gray(image), se);
        BW = imbinarize(im,'global');
        BW = bwareaopen(BW, aoi.radius, 8);
        stats = regionprops(BW,im,'WeightedCentroid','MaxIntensity', 'Area');
end

centroids = cat(1,stats.WeightedCentroid);
numROIs = size(centroids,1);

% Remove spots based on distance
spotsRemoved =[];
for i = 1:numROIs
    dist = sqrt((centroids(i,1) - centroids(:,1)).^2 + (centroids(i,2) - centroids(:,2)).^2);
    removeIndex = find(dist < aoi.minSeperation);
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
boundingBox = makeBoundingBox(centroids, aoi.area);

if aoi.gauss
    aoiSigma = centroids*0;
    gauss = cell(numROIs,1);
    wbGauss = waitbar(0, 'Gaussian refinement...');
    for i = 1:numROIs
        c = centroids(i,:) - boundingBox(i,1:2)+1;
        imageTemp = imcrop(image, boundingBox(i,:));
        
        params0 = [c(1), c(2), aoi.radius-1, aoi.radius-1, imageTemp(round(c(1)), round(c(2))), 0];
        [params,~] = fitgaussian2d(imageTemp, params0, 2, aoi.tol);
        newCentroid = params(1:2);
        dist = calcDist(c, newCentroid);
        if dist <= aoi.radius
            centroids(i,1)= params(1) + boundingBox(i,1)-1;
            centroids(i,2)= params(2) + boundingBox(i,2)-1;
            boundingBox(i,:) = makeBoundingBox(centroids(i,:), aoi.area);
            aoiSigma(i,:) = params(3:4);
        else
            % Remove the AOI? Probably empty or multiple occupied
        end
        waitbar(i/numROIs, wbGauss);
    end
    close(wbGauss)
end

% update the bounding box
boundingBox = makeBoundingBox(centroids, aoi.area);

% Make roi struct
rois = struct;
for i = 1:numROIs
    rois(i,1).Centroid = centroids(i,:);
    % rois(i,1).sigma = aoiSigma(i,1); % they are the same value
    rois(i,1).boundingBox = boundingBox(i,:);
    rois(i,1).pixelList = boundBoxToPixels(rois(i,1).boundingBox );
end

if showResult
    autoImageBC(image, 1); hold on
    for i = 1:length(rois)
        rectangle('Position', rois(i).boundingBox(1,:),...
            'EdgeColor','r', 'LineWidth', 1, 'Curvature',[0 0]);
    end
    scatter(centroids(:,1), centroids(:,2), 'y*')
    pause(1)
end

