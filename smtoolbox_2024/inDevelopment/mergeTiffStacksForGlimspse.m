%% Merge 2 tiff stacks into single image for imscroll. 
% David S. White 
% 2021-11-18 
% MIT

%% Load the left (633 nm) image stack. 
[fileLeft, pathLeft] = uigetfile('.tif');
videoLeft = loadtiff([pathLeft, pathRight]); 

%% Load the left (532 nm) image stack
[fileRight, pathRigth] = uigetfile('.tif');
videoRight = loadtiff([pathRigth, fileRight]); 

%% Should be the same size (no error checking). Stack horizontally 
videoStacked = [videoLeft, videoRight]; 

%% Now save the new image stack as tif
[saveFile, savePath] = uiputfile('.tif'); 
saveastiff(videoStacked, [savePath, saveFile]); 

