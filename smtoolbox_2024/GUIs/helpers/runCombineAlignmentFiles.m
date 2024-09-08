%% run combineAlignmentFiles.m
% David S. White
% 2023-02-07

% select files
[files,paths] = getMultipleFiles();
% need debug if only one file selected... but why would one do that? 
nFiles = length(files); 

alignTransformIn = cell(nFiles,1);
for i = 1:nFiles
    load([paths{i},files{i}]);
    alignTransformIn{i} = alignTransform;
end

% merge all
alignTransform = combineAlignmentFiles(alignTransformIn);
