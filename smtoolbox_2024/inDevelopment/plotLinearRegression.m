function [h, plot_data] = plotLinearRegression(mdl, x, y, varargin)
% David S. White 
% 2023-09-08

% Function to plot results of linearRegression.m
% mdl = result from linearRegression 

% Defaults
parent = [];
showConfidenceInterval = true; 
confidenceInterval = 0.95;
color = 'r';
marker = 'o';
markerSize = 5; 
markerFaceColor = 'w'; 
markerEdgeColor = 'k';
xmin = min(x);
xmax = max(x);

for i = 1:2:length(varargin)-1
    switch varargin{i}
        case 'parent'
            % figure handle
            parent = varargin{i+1};
        case 'showConfidenceInterval'
            showConfidenceInterval = varargin{i+1};
        case 'confidenceInterval'
            confidenceInterval = varargin{i+1};
        case 'color'
            color = varargin{i+1};
        case 'marker'
            marker = varargin{i+1};
        case 'markerSize'
            markerSize = varargin{i+1};
        case 'markerFaceColor'
            markerFaceColor = varargin{i+1};
        case 'markerEdgeColor'
            markerEdgeColor = varargin{i+1};
        case 'xmin'
            xmin = varargin{i+1};
        case 'xmax'
            xmax = varargin{i+1};
    end
end

% adjust data
x = x(:);
y = y(:);

% evaluate over xmin:xmax
xx = linspace(xmin, xmax, 100);
xx = xx(:);
yy = feval(mdl, xx);

% 95% Confidence intervals
if strcmp(get(parent,'type'),'figure')
    h = parent;
else
    h = figure;
end
hold on;
if showConfidenceInterval
    CIF = predint(mdl, xx, confidenceInterval, 'Functional');
    patch([xx' fliplr(xx')], [CIF(:,2)' fliplr(CIF(:,1)')],...
        color, 'FaceColor', color, 'EdgeColor','none', 'FaceAlpha', 0.1);
end
plot(xx, yy, '-', 'color', color)
plot(x, y, 'Marker', marker, 'MarkerSize', markerSize, 'MarkerFaceColor', markerFaceColor, 'MarkerEdgeColor', markerEdgeColor, 'Linestyle', 'None');
xlabel('x')
ylabel('y')
xlim([xmin, xmax])

plot_data.x = x; 
plot_data.y = y; 
plot_data.xx = xx; 
plot_data.yy  = yy; 
plot_data.CI = CIF;

