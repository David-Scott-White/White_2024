function [videos, time_s] = selectBioformat(videos, time_s)
% David S. White
% 2021-09-14

% Pop-up to select which videos to retain (assumes two channels); 

% Launch a quick GUIDE app (should make programically...) 


h = guiSelectChannel()


