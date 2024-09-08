function [h,sortedIdx] = plotRaster(events, varargin)
% -------------------------------------------------------------------------
% Raster plot of idealized traces
% -------------------------------------------------------------------------
%
% events = cell(N,1) events for each molecule. see findEvents.m
%
% Example Code (Temp)
%
% N = 100; 
% events = cell(N,1);
% for i = 1:N
%     fps = 1;
%     q = [0, 0.0025; 0.05, 0];
%     ss = simulateQ(q, [0, 1], 1e2, fps);
%     events{i} = findEvents(ss); 
% end
% sortTraces = 1; 
% plotRaster(traces, 'frameRate_s', 1/fps)
%
% Author: David S. White
% Updated: 2021-10-28
% License: MIT
% -------------------------------------------------------------------------
%%
% Set Default parameters
sortTraces  = 1;   % boolean (ascending for now) 
colorScheme = [0.7,0.7,0.7]; % need a gradient for multiple states
frameRate_s = 1;   % else, returns frames
modFrameRate = 0;
numMolecules = length(events); 
plotFractionBound = 1; 
fractionBoundColor = [0.8500, 0.3250, 0.0980]; 
seed = 42; % for random number generator in numEvents
width = 1; 
for i = 1:2:length(varargin)-1
    switch varargin{i}
        case 'sortTraces'
            sortTraces = varargin{i+1};
        case 'colorScheme'
            colorScheme = varargin{i+1};
        case 'frameRate_s'
            frameRate_s = varargin{i+1};
            modFrameRate = 1;
        case 'plotFractionBound'
            plotFractionBound = varargin{i+1};
        case 'fractionBoundColor'
            fractionBoundColor = varargin{i+1};
        case {'numMolecules', 'N'}
            % randomly show N number of traces of the full set
            numMolecules = varargin{i+1};
        case 'seed'
            seed = varargin{i+1};
        case 'width'
            width = varargin{i+1};
    end
end

concatEvents = cell2mat(events);
N0 = length(events);
minState = min(concatEvents(:,4));
maxFrame = max(concatEvents(:,2));

% Compute fraction bound
if plotFractionBound
    fractionBound = zeros(maxFrame,1);
    for i = 1:length(events)
        if ~isempty(events{i})
            sequence = eventsToSequence(events{i});
            fractionBound = fractionBound + sequence;
        end
    end
    fractionBound = fractionBound/length(events);
else
    fractionBound = [];
end

% use all events or subsample ? 
if numMolecules < length(events)
    rng(seed, 'twister');
    idx1 = randperm(length(events));
    idx2 = idx1(1:numMolecules);
    events = events(idx2); 
else
    numMolecules = length(events);
end


% Sort traces?
N = size(events,1);
if sortTraces
    firstEventFrame = zeros(N,1);
    for i = 1:N
        if isempty(events{i})
            firstEventFrame(i) = maxFrame+1;
        elseif events{i}(1)==0
            firstEventFrame(i) = maxFrame+1;
        elseif size(events{i},1) == 1
            if events{i}(4) == minState
                firstEventFrame(i) = events{i}(2);
            else
                firstEventFrame(i) = events{i}(1);
            end
        else
            if events{i}(1,4) > min(events{i}(:,4))
                firstEventFrame(i) = 1;
            else
                firstEventFrame(i) = events{i}(2,1);
            end
        end
    end
    [~,sortedIdx] = sort(firstEventFrame);
    events = events(sortedIdx, :);
end



% build figure
h = figure;
ylabel('Molecule')
hold on
xlim([0, maxFrame*frameRate_s]);
ylim([1, N])
for i = 1:N
    if ~isempty(events{1})
        for j = 1:size(events{i},1)
            if events{i}(j,4) > minState
                s0 = events{i}(j,1)*frameRate_s;
                s1 = events{i}(j,3)*frameRate_s;
                % xstart, ystart, xlength, ylength
                rectangle('Position',[s0, i-0.5, s1, width], 'EdgeColor', colorScheme, 'FaceColor', colorScheme);
            end
        end
    end
end
h.Children(1).YColor = 'k';

% should add text to top left for number of molecuels in plot 
h.Children(1).YColor = 'k';
if ~isempty(fractionBound)
         yyaxis right
         ylabel('Fraction Bound')
    hold on
    if modFrameRate
        time_s = frameRate_s:frameRate_s:(length(fractionBound)*frameRate_s);
        plot(time_s, fractionBound, '-', 'color', fractionBoundColor, 'linewidth',1);
    else
        plot(fractionBound, '-', 'color', fractionBoundColor, 'linewidth',1);
    end
    if sum(fractionBound>1)
        ylim([-1, N])
        ylabel('Molecules Bound')
    else
        ylim([0,1])
        ylabel('Fraction Bound')
    end
    h.Children(1).YColor = fractionBoundColor;
end

if modFrameRate
    xlabel('Time (s)');
else
    xlabel('Frames');
end
set(gca, 'tickdir', 'out');

xL=xlim;
yL=ylim;
text(xL(1),0.99*yL(2),['N = ', num2str(N0)],'HorizontalAlignment','left','VerticalAlignment','top', 'FontName', 'Arial')

