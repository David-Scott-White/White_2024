function exportInfo = getExportVideoInfo()
% Get values for writeImagesToVideo.m
% 
% David S. White
% 2022-02-01
% MIT
exportInfo = struct; 
prompt = {'Channel', 'FrameRate (s)', 'Show AOI (bool)',...
    'Show Time (bool)', 'Scale Bar (um)', 'Pixels per um', 'Start Frame', 'Stop Frame'};
dlgtitle = 'Export Video Input';
dims = [1 35];
definput = {'1', '10', '1', '1', '10', '0.216', '1', ''};
answer = inputdlg(prompt,dlgtitle,dims,definput);

% emperical
if ~isempty(answer)
    exportInfo.channel = str2double(answer{1});
    exportInfo.frameRate_s = str2double(answer{2});
    exportInfo.showAOI = str2double(answer{3});
    exportInfo.showTime = str2double(answer{4});
    exportInfo.scaleBar_um = str2double(answer{5});
    exportInfo.umPerPixel = str2double(answer{6});
    exportInfo.startFrame = str2double(answer{7});
    exportInfo.stopFrame = str2double(answer{8});
end
end