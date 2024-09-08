function guiSaveData(allData)
% get desired variables
data = struct;
if isfield(allData,'rois')
    data.rois = allData.rois;
end
if isfield(allData,'time')
    data.time = allData.time;
end
if isfield(allData,'videoPath')
    data.videoPath = allData.videoPath;
end
if isfield(allData,'aoiInfo')
    data.aoiInfo = allData.aoiInfo;
end
if isfield(allData,'driftList')
    data.driftList = allData.driftList;
end
if isfield(allData,'alignTransform')
    data.alignTransform = allData.alignTransform;
end
if isfield(allData,'numPixels')
    data.numPixels = allData.numPixels;
end
if isfield(allData,'aoiInfo')
    data.aoiInfo = allData.aoiInfo;
end

tempName = [allData.videoPath, 'data.mat'];
[file, path] = uiputfile(tempName);
disp('> Saving data...')
save([path,file], 'data', '-v7.3');
msgbox(' > Data saved:');
disp(['  >> Saved:', [path, file]]);
