function exportFrame(image, varargin)
% -------------------------------------------------------------------------
% Export frame from video in smViewer to publishFigure (.mat, .png, .pdf)
%
%
% David S. White
% 2021-02-01
% MIT
% -------------------------------------------------------------------------
aoi = 0; % bounding box only
scaleBar = 0;
umPerPixel = 0;
time_s = 0;
filePath = 0;
scaleLimits =[];
fontSize = 18;

% add option for figure size in future, display range

for i = 1:2:length(varargin)-1
    switch varargin{i}
        case 'aoi'
            aoi = varargin{i+1};
        case 'scaleBar'
            scaleBar = varargin{i+1};
        case 'umPerPixel'
            umPerPixel = varargin{i+1};
        case 'format'
            format = varargin{i+1};
        case 'quality'
            quality = varargin{i+1};
        case 'time_s'
            time_s = varargin{i+1};
        case 'filePath'
            filePath = varargin{i+1};
        case 'scaleLimits'
            scaleLimits = varargin{i+1};
    end
end

if ~filePath
    [file, path] = uiputfile('video.mp4');
    filePath = [path,file];
end

%%
figExist = findobj('type','figure','name','Frame Export');
if ~isempty(figExist)
    close('Frame Export');
end

[xImage, yImage] = size(image);
if isempty(scaleLimits)
    [x1, x2] = autoImageBC(image);
else
    x1 = scaleLimits(1); 
    x2 = scaleLimits(2);
end

hFig = figure('Name', 'Frame Export', 'NumberTitle', 'off', 'Visible', 'off');
imshow(image, 'DisplayRange', [x1,x2],'Border', 'tight');

hFig.Position(1) = 0;
hFig.Position(2) = 0;
hFig.Position(3) = xImage;
hFig.Position(4) = yImage;
hFig.Visible = 'on';
        
pause(1)
hold on;
if aoi
    for j = 1:length(aoi)
        rectangle('Position', aoi(j,:), 'EdgeColor', 'r', 'linewidth', 0.5);
    end
    pause(1);
end

if time_s
    s = seconds(round(time_s));
    s.Format = 'hh:mm:ss';
    timestr = char(s);
    text(10,10, timestr,'Color','w','FontSize', fontSize, 'FontName', 'Helvetica');
    pause(1)
end

if scaleBar && umPerPixel
    pixels = scaleBar/umPerPixel;
    p0 = xImage-10;
    p1 = xImage-10-pixels;
    sbX = p1:p0;
    sbY = sbX*0 + yImage-10;
    plot(sbX, sbY,'-w','linewidth', 3)
    text(p0, yImage-25, [num2str(scaleBar), '\mum'], 'Color', 'w','FontSize', fontSize, 'HorizontalAlignment', 'right','FontName','Helvetica');
end

[file,path] = uiputfile([filePath, '.pdf']);
if file
    filePathSave = [path,file(1:end-4)];
    % hFig.Children.Children = flip(hFig.Children.Children);
    publishFigure(hFig, 'figureName', filePathSave, 'saveOnly', 1);
    
    % note, text will be behind image in PDF (do not know how to fix...)
    % and wrong color... and wrong font. WTF... 
    
end
close('Frame Export');