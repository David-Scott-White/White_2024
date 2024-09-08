function AOIs = transformAOI(AOIs, radius, channelTransformation)
% David S. White 
% 2021-09-19 

% Transform AOIs using a provided tform
%
% rois: struct with rois.centroid
% tform: result from image regresitation 
% radius 

AOIs = [AOIs, AOIs];
if ~exist('channelTransformation', 'var') || isempty(channelTransformation)
     msgbox('Warning: No Transform provided.');
else
    if isfield(AOIs, 'Centroid')
        centroids1 = vertcat(AOIs(:,1).Centroid);
    elseif isfield(AOIs, 'centroid')
        centroids1 = vertcat(AOIs(:,1).centroid);
    end
    % centroids = applyTransform(centroids, channelTransformation);
    centroids2 = transformPointsInverse(channelTransformation.tformAOI, centroids1);
    diameter = AOIs(1,1).boundingBox(4); 
    boundingBox = makeBoundingBox(centroids2, diameter); 
    for i = 1:length(AOIs)
        AOIs(i,2).centroid = centroids2(i,:);
        AOIs(i,2).Centroid = centroids2(i,:);
        AOIs(i,2).boundingBox = boundingBox(i,:);
    end
end

