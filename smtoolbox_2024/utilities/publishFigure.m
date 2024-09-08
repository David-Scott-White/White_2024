function h = publishFigure(h, varargin)
% -------------------------------------------------------------------------
% Modify Figure (h) for publication
% -------------------------------------------------------------------------
%
%
% David S. White
% 2021-10-19
% -------------------------------------------------------------------------

% WORK IN  PROGRESS
% set defaults

saveOnly = 0;
figureName = [];
figureSize = [];
fontName = 'Arial';
fontSize = 7;
legendMarkerSize = 5;
legendFontSize = 5;
legendTokenSize = [10, 6];
tickDir = 'out';
box = 'off';
plotColor = 'k';
showImage = 0;
tickLength = 0.02;
fileType = {'.pdf', '.png'};

for i = 1:2:length(varargin)-1
    switch varargin{i}
        case 'figureName'
            figureName = varargin{i+1};
        case 'figureSize'
            figureSize = varargin{i+1};
        case 'fontName'
            fontName = varargin{i+1};
        case 'fontSize'
            fontSize = varargin{i+1};
        case 'tickDir'
            tickDir = varargin{i+1};
        case 'box'
            box = varargin{i+1};
        case 'plotColor'
            plotColor = varargin{i+1};
        case 'showImage'
            showImage = varargin{i+1};
        case 'legendMarkerSize'
            legendMarkerSize = varargin{i+1};
        case 'legendFontSize'
            legendFontSize = varargin{i+1};
        case 'legendTokenSize'
            legendTokenSize = varargin{i+1};
        case 'saveOnly'
            saveOnly = varargin{i+1};
        case 'tickLength'
            tickLength = varargin{i+1};
        case 'fileType'
            fileType = varargin{i+1};
    end
end


if ~saveOnly
    
    set(gca,...
        'fontsize',fontSize,...
        'fontname', fontName,...
        'tickdir',tickDir,...
        'box',box,...
        'TickLength', [tickLength, tickLength],...
        'XColor', plotColor);
    
    set(h, 'DefaultAxesFontName', fontName);
    set(h, 'DefaultTextFontName', fontName);
    
    if ~isempty(figureSize)
        h.Units = 'inches';
        h.Position = [0, 0, figureSize(1), figureSize(2)];
    end
    
    hh  = findobj('-property','FontSize');
    for i = 1:length(hh)
        hh(i).FontSize = fontSize;
    end
    
    % check for legend, modify
    if ~isempty(legendTokenSize)
        hh  = findobj('Type', 'Legend');
        for i = 1:length(hh)
            hh(i).FontSize = legendFontSize;
            hh(i).ItemTokenSize = legendTokenSize;
        end
    end
end

% move axes on top of data
set(gca, 'Layer', 'top')

% Save (as a .mat, .pdf, and .png);
if isempty(figureName)
    [file,path] = uiputfile('.fig');
    figureName = [path,file(1:end-4)];
end
% savefig(h, figureName);
% print(h, figureName,'-dpdf',  '-r0')
% OMG THIS FINALLY WORKS. Even Exports Arial Font Now!. Just a bit slow
if any(strcmp(fileType,'.pdf'))
    exportgraphics(h,[figureName, '.pdf'], 'ContentType','vector');
elseif any(strcmp(fileType,'.png'))
    exportgraphics(h, [figureName, '.png'], 'Resolution', 300);
elseif any(strcmp(fileType,'.fig'))
    savefig(h, figureName);
end

% Note may need to increase java heap memory. Prefences->General->Java Heap
% Memory

if showImage
    close(h);
    A = imshow(imread([figureName,'.png']), 'InitialMagnification', 100);
end

end
