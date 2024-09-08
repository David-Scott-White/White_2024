%% Launch smViewer
% -------------------------------------------------------------------------
% Launch the smViewer software. 
% Loads the Home and creates "smData" variable to be shared across GUI. 
%
% To do: 
% > Create a class for data (or rather, all the var experiment, aoi, etc..)
%
% David S. White 
% 2021-09-02
% -------------------------------------------------------------------------
% intit all fields here
data = struct;
% data.vp = struct; % parameters for video processing
% data.tv = struct; % parameters for trace viewer 

smHome(data)
clear data