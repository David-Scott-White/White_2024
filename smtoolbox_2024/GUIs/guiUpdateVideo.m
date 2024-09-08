function handles = guiUpdateVideo(hObject, handles, videoIndex, updateAOI, updateImage)
% -------------------------------------------------------------------------
% Update the video image shown in smVideoProcessing.
% -------------------------------------------------------------------------
%
% run  after click of sliders and loading videos
%
% David S. White
% dwhite7@wisc.edu
% GNU General Public License v3 (2021-09-27)
% -------------------------------------------------------------------------
if nargin < 4
    updateAOI = 0;
end
if nargin < 5
    updateImage = 0;
end

if isempty(handles.data.videos{videoIndex})
    return
end

frame = handles.data.frame{videoIndex}(1);
image = handles.data.videos{videoIndex}(:,:,frame);
switch videoIndex
    case 1
        figureHandle = handles.axes1;
        displayRange = [handles.slider_BCmin1.Value, handles.slider_BCmax1.Value];
        handles.slider_BCmin1.Max = handles.slider_BCmax1.Value-1;
        handles.slider_BCmax1.Min = handles.slider_BCmin1.Value+1;
        showBoundingBox = handles.showAOI1.Value;
        handles.frame1.String = num2str(frame);
        handles.text_bcmin1.String = num2str(round(displayRange(1)));
        handles.text_bcmax1.String = num2str(round(displayRange(2)));
    case 2
        figureHandle = handles.axes2;
        displayRange = [handles.slider_BCmin2.Value, handles.slider_BCmax2.Value];
        handles.slider_BCmin2.Max = handles.slider_BCmax2.Value-1;
        handles.slider_BCmax2.Min = handles.slider_BCmin2.Value+1;
        showBoundingBox = handles.showAOI2.Value;
        handles.frame2.String = num2str(frame);
        handles.text_bcmin2.String = num2str(round(displayRange(1)));
        handles.text_bcmax2.String = num2str(round(displayRange(2)));
end

% Display the image
showBB = 0;
if updateImage
    imshow(image,'DisplayRange', displayRange, 'Parent', figureHandle);
else
    if isempty(figureHandle.Children)
        imshow(image,'DisplayRange', displayRange, 'Parent', figureHandle);
        if isfield(handles.data, 'rois')
            if ~isempty(handles.data.rois) && showBoundingBox
                showBB = 1;
            end
        end
    else
        figureHandle.Children(end).CData = image;
        if ~isempty(handles.data.rois) && showBoundingBox
            showBB = 1;
        end
    end
end
if showBB
    if length(figureHandle.Children) == 1
        for i = 1:size(handles.data.rois,1)
            bb = handles.data.rois(i,videoIndex).boundingBox(frame,:);
            % change appearance to enclose pixels, rather than cut through
            % them
            bb(1) =bb(1)-0.5; 
            bb(2) = bb(2)-0.5;
            bb(3:4) = bb(3:4)+1;
            rectangle('Parent',figureHandle, 'Position', bb,...
                'EdgeColor','r', 'LineWidth', 1);
        end
    elseif updateAOI
        delete(findobj(figureHandle.Children, 'type', 'Rectangle'));
        for i = 1:size(handles.data.rois,1)
            bb = handles.data.rois(i,videoIndex).boundingBox(frame,:);
            % change appearance to enclose pixels, rather than cut through
            % them
            bb(1) =bb(1)-0.5; 
            bb(2) = bb(2)-0.5;
            bb(3:4) = bb(3:4)+1;
            rectangle('Parent',figureHandle, 'Position',bb,...
                'EdgeColor','r', 'LineWidth', 1);
        end
    end
else
    delete(findobj(figureHandle.Children, 'type', 'Rectangle'));
end
drawnow();

end