function [files, paths] = getMultipleFiles(ext)
% -------------------------------------------------------------------------
% Get File Paths of Multiple Files Across Different Directories
% -------------------------------------------------------------------------
% Continuous selection of different files across different directories.
% Useful replacement for "Multiselect" in uigetfile.
%
% Input:
%   N/A
%
% Output:
%   files = cell(N); selected files
%   paths = cell(N); selected paths
%
% Author: David S. White
% Upated: 2021-10-08
% License: MIT
% -------------------------------------------------------------------------
files = {};
paths = {};
nFiles = 1;
if nargin < 1
    ext = '*.mat';
end
[files{nFiles}, paths{nFiles}] = uigetfile(ext, 'Select a File');
answer = questdlg('Select another file?', 'User Answer', 'Yes', 'No', 'No');
disp('Files selected: ')
disp(['>> ', files{1}])
while strcmp(answer, 'Yes')
    nFiles = nFiles + 1;
    [files{nFiles}, paths{nFiles}] = uigetfile(ext, 'Select a File');
    disp(['>> ', files{nFiles}])
    answer = questdlg('Select another file?', 'User Answer', 'Yes', 'No', 'No');
end



