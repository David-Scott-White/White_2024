function Intervals = concatIntervals(sameDirectory)
% -------------------------------------------------------------------------
% Combine Intervals.mat files into a single file (imscroll)
% -------------------------------------------------------------------------
% Input:
%   sameDirectory. boolean.
%
% Output:
%   Combined "Intervals" structure read by imscroll. Saves the file.
%
% Requires:
%   getMultipleFiles.m (DSW)
%
% Author: David S. White
% Updated: 2021-10-08
% License: MIT
% -------------------------------------------------------------------------
if ~exist('sameDirectory','var') || isempty(sameDirectory)
    sameDirectory = 0;
end
if sameDirectory
    [files, paths] = uigetfile('*.mat', 'MultiSelect','on');
else
    [files, paths] = getMultipleFiles();
end
nFiles = length(files);
for i = 1:nFiles
    if sameDirectory
        temp = load([paths, files{i}]);
    else
        temp = load([paths{i}, files{i}]);
    end
    if i == 1
        Intervals = temp.Intervals;
    else
        Intervals.AllTracesCellArray = [Intervals.AllTracesCellArray; temp.Intervals.AllTracesCellArray];
        Intervals.CumulativeIntervalArray = [Intervals.CumulativeIntervalArray; temp.Intervals.CumulativeIntervalArray];
    end
end
[saveFile, savePath] = uiputfile('Intervals.mat');
if saveFile
    save([savePath, saveFile], 'Intervals')
end
end


