function guiUpdateAlignImage(handles)
% David S. White 
% 2021-09-19

% Quick function to update align video

video1 = handles.data.alignVideos{1}(:,:,1); % adjust for mean value...
video2 = handles.data.alignVideos{2}(:,:,1);
[x1, x2] = autoImageBC(video1); 
[y1, y2] = autoImageBC(video2); 

imshow(video1, 'DisplayRange', [x1,x2], 'Parent', handles.axes_channel1);
imshow(video2, 'DisplayRange', [y1,y2], 'Parent', handles.axes_channel2);

imshowpair(imadjust(mat2gray(video1)), imadjust(mat2gray(video2)), 'Parent', handles.axes_result);  