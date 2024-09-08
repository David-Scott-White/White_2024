function errorPlot(x, y, yE, varargin)
% -------------------------------------------------------------------------
% Shaded Errorbar Plot
% -------------------------------------------------------------------------
% David S. White
% dwhite7@wisc.edu
% 2021-02-07
%
% Inputs:
% -------
% x  = [1xN] array of X Values
% y  = [1xN] array of X Values
% yE = [1xN] array of error for Y. Currently does not allow for different
%      lower and upper limits.
%
% Output:
% -------
%
%
% Notes:
% -------
%
% -------------------------------------------------------------------------
%` Check input, set defaults if not pro
[~,n] = size(x); if n == 1; x = x'; end
[~,n] = size(y); if n == 1; y = y'; end
[m,n] = size(yE); if m > n; yE = yE'; end
[m,~] = size(yE);
if m == 1
    yL = y-yE; % lower bound
    yU = y+yE; % upper bound
elseif m == 2
    yL = yE(1,:); 
    yU = yE(2,:);
end

%% Set plot parameters (default values if not provided)
plotColor = [];
lineWidth = 0.5;
lineColor = 'r';
edgeColor = 'r';
patchColor = 'r';
patchAlpha = 0.05;
mainMarker = ':';
edgeMarker = '-';
mainFaceColor = 'w';
mainEdgeColor = 'k';
mainMarkerSize = 5; 

% check varargin and change plot parameters.
for i = 1:2:length(varargin)-1
    switch varargin{i}
        case 'lineWidth'
            lineWidth = varargin{i+1};
        case 'lineColor'
            lineColor = varargin{i+1};
        case 'edgeColor'
            edgeColor = varargin{i+1};
        case 'mainMarker'
            mainMarker = varargin{i+1};
        case 'edgeMarker'
            edgeMarker = varargin{i+1};
        case 'patchColor'
            patchColor = varargin{i+1};
        case 'patchAlpha'
            patchAlpha = varargin{i+1};
        case 'mainMarkerSize'
            mainMarkerSize = varargin{i+1};
        case 'plotColor'
            plotColor = varargin{i+1};
    end
end
if ~isempty(plotColor)
    lineColor = plotColor;
    edgeColor = plotColor;
    patchColor = plotColor;
end

%% Plot
patch([x fliplr(x)], [yU fliplr(yL)], edgeMarker, 'FaceColor', patchColor, 'EdgeColor','w', 'FaceAlpha', patchAlpha); hold on
plot(x, y,  mainMarker, 'color', lineColor,'MarkerFaceColor', mainFaceColor,...
    'MarkerEdgeColor',mainEdgeColor, 'MarkerSize', mainMarkerSize, 'linewidth',lineWidth);
plot(x, yL, edgeMarker, 'color', edgeColor, 'linewidth', lineWidth);
plot(x, yU, edgeMarker, 'color', edgeColor, 'linewidth', lineWidth);

end