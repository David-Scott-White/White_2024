function alignTransform = combineAlignmentFiles(alignTransformIn)
%% Create an ensemble alignment file from multiple alignment files
% David S. White 
% 2022-12-05
%
% Patch onto smViewer. Have to call outside

%  Input.
% alignTransformIn cell array of alignTransform generated in smVideoProcessing

% Output. 
% alignTransform that considers all datapoints from all
% alignTransformIn

% assume similarity transform is fine for each
alignTransform = alignTransformIn{1};

n = length(alignTransformIn);
for i = 2:n % first one already copied
    alignTransform.tformImage.T(:,1:2) = alignTransform.tformImage.T(:,1:2) + alignTransformIn{i}.tformImage.T(:,1:2);
    
    alignTransform.centroids{1} = [alignTransform.centroids{1}; alignTransformIn{i}.centroids{1}];
    alignTransform.centroids{2} = [alignTransform.centroids{2}; alignTransformIn{i}.centroids{2}];
end

% average image alignment (similarity)
alignTransform.tformImage.T(:,1:2) = alignTransform.tformImage.T(:,1:2)/n;

% new aoi alignment (affine)
nCentroids = length(alignTransform.centroids{1});
alignTransform.tformAOI = fitgeotrans(alignTransform.centroids{2}, alignTransform.centroids{1}, 'affine');
centroids3 = transformPointsInverse(alignTransform.tformAOI, alignTransform.centroids{1});

alignTransform.rmse = sqrt(sum((centroids3-alignTransform.centroids{2}).^2)./nCentroids);
disp(['RMSE (x, y): ', num2str(alignTransform.rmse)]);

[saveFile,savePath] = uiputfile('combinedAlignmentFile.mat')
save([savePath,saveFile], 'alignTransform');



