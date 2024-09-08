function [alignParm, image3] = batchAlignImagesAOI(image1, image2, aoi, flip2)
% Align Two detectors using similarity & affine transform 
% Overview: 
%   For each pair of images, compute a similarity transform. Store. 
%   Auto detect spots in each image following input from aoi (findAOI.m) 
%   Auto match spots into pairs by min dist, remove outliers & store
%   For all pairs of spots across all images, compute an affine transform. 
% 
%
% Input: 
%   image1 = cell of [x,y,1] data; 
%   image2 = cell of [x,y,1] data
%   aoi = aoi parms. see findAOI.m
%   flip2 = bool (& temp var). flip image 1 horizonataly.         
%
% Output
%   alignParm = struct; 
%       alignParm.tformImage
%       alignParm.tformSpots
%       alignParm.rmse
%
%
% Requires Image Processing Toolbox
%
% David S. White 
% 2022-02-03
% License: MIT
% -------------------------------------------------------------------------

if nargin < 3 || isempty(aoi)
    aoi = struct;
    aoi.area = 3;
    aoi.radius = 2;
    aoi.fp = 20;
    aoi.gauss = 1;
    aoi.tol = 1e-5;
    aoi.method = 'GLRT';
    aoi.minSeperation = aoi.radius;
end

if nargin < 4
    flip2 = 0;
end

if isempty(image1)
    [files, paths] = getMultipleFiles('*.tif');
    nFiles = length(files); 
    image1 = cell(nFiles,1); 
    image2 = cell(nFiles,1); 
    for i = 1:nFiles
        videos = loadBioFormat(files{i}, paths{i}, 1, 0);
        image1{i} = mean(videos{end},3);
        image2{i} = mean(videos{1},3);
        if flip2
            image2{i} = flip(image2{i},2);
        end
    end
else
    nFiles = length(image1);
end

%% Compute similarity transform between all images. Use average
refObj1 = imref2d(size(image1{1}));
if size(image1{1},3) == size(image2{1}, 3)
    refObj2 = refObj1;
else
    refObj2 = imref2d(size(image2{1}));
end
method = 'zscore';
allImageTform = zeros(3,3,nFiles);
for i = 1:nFiles
    image1s = im2single(scaleFeature(image1{i}, method));
    image2s = im2single(scaleFeature(image2{i}, method));
    
    tform = imregcorr(image2s,refObj2,image1s,refObj1,...
        'transformType', 'similarity', 'Window',false);
    allImageTform(:,:, i) = tform.T;
    
end
imageTform = tform;
imageTform.T = median(allImageTform,3); % mean or median transform? 
image3 = cell(nFiles,1);
for i = 1:nFiles
    image3{i} = imwarp(image2{i}, refObj2, imageTform, 'OutputView', refObj1, 'SmoothEdges', true);
end

%% Auto detect spots in each image1 & image3.  (no longer using image 2)
centroids1 = cell(nFiles,1); 
centroids2 = cell(nFiles,1);
for i = 1:nFiles
    [centroids1{i}, centroids2{i}] = localGetPairs(image1{i}, image3{i}, aoi); % need to replace with aoi
end

%% fit all centroids with an affine transform
c1 = cell2mat(centroids1);
c2 = cell2mat(centroids2);

nCentroids = length(c1); 
disp(['Number of Centroid Pairs:', num2str(nCentroids)]);

tform2 = fitgeotrans(c2, c1, 'affine');
c3 = transformPointsInverse(tform2, c1);

rmse = sqrt(sum((c3-c2).^2)./nCentroids);
disp(['RMSE (x, y): ', num2str(rmse)]);
% disp(['Error (nm): ', num2str(rmse2*260)]); % have an info output when
% building experiments/set ups is an option in the full software

% output 
alignParm.tformImage = imageTform;
alignParm.tformSpots = tform2; 
alignParm.rmse = rmse;

end

%% 
function [centroids1, centroids2] = localGetPairs(im1, im2, aoi)
% Local function to match pairs between the images
if isempty(aoi)
    aoi = struct;
    aoi.area = 3;
    aoi.radius = 2;
    aoi.fp = 20;
    aoi.gauss = 1;
    aoi.tol = 1e-5;
    aoi.method = 'GLRT';
    aoi.minSeperation = aoi.radius;
end

rois1 = findAOI(im1, aoi, 0);
rois2 = findAOI(im2, aoi, 0);
centroids1 = vertcat(rois1.Centroid);
centroids2 = vertcat(rois2.Centroid);
pairs = centroids1*0;
keep = ones(length(centroids1),1);
for i = 1:length(centroids1)
    dist = sqrt((centroids1(i,1) - centroids2(:,1)).^2 + (centroids1(i,2) - centroids2(:,2)).^2);
    [minDist, j] = min(dist);
    if minDist > 2
        keep(i) = 0;
    end
    if sum(pairs(:,2)==j)
        keep(i) = 0;
    else
        pairs(i,1) = i;
        pairs(i,2) = j;
    end
end
pairs = pairs(keep==1, :);

centroids1 = centroids1(pairs(:,1),:);
centroids2 = centroids2(pairs(:,2),:);
dist = sqrt((centroids1(:,1) - centroids2(:,1)).^2 + (centroids1(:,2) - centroids2(:,2)).^2);
outlierPresent = 1;
while outlierPresent
    idx = find(isoutlier(dist)==1);
    if ~isempty(idx)
        centroids1(idx,:) = [];
        centroids2(idx,:) = [];
        dist = sqrt((centroids1(:,1) - centroids2(:,1)).^2 + (centroids1(:,2) - centroids2(:,2)).^2);
    else
        outlierPresent = 0;
    end
end
end
    