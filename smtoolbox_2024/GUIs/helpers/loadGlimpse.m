function [path, video, time_s] = loadGlimpse()

path = uigetdir();
splitImage = 0;
[images, time_s] = glimpseToTif(path, splitImage);
video = cell(1,1);
video{1} = images;