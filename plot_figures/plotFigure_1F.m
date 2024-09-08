%% Plot Figure 01-v2
% David S. White
% 2024-01-01

% add White_2023 to file path
close all 
clear

fontSize = 7; 
figureSize = [5,1]; 
figurePath = '/Users/work/Documents/Research/White_2023/figures-raw/figure01/';

%% Store output 
data_out = {}

%% Example Trace of WT and -1A 
load('20230616_img01_9bp-1A-1nM_U1C-100nM_DMSO.mat')


% plot example trace 
i = 648; 

h2 = figure; hold on
yline(data.rois(i,2).fit.components(1,2), '--', 'color', [0.5, 0.5, 0.5])
yline(data.rois(i,2).fit.components(2,2), '--', 'color', [0.5, 0.5, 0.5])
plot(data.time{2}(:,end), data.rois(i,2).timeSeries, 'color', [0.4660 0.6740 0.1880])
plot(data.time{2}(:,end), data.rois(i,2).fit.ideal, '-k')
xlim([0, 1800]);
xticks(0:300:1800);
xlabel('Time (s)');
ylim([0, 700])
yticks(0:200:800)
% ylabel('Fluorescence (AU)')
publishFigure(h2,...
    'figureName', [figurePath, 'trace_9bp-1A_1nM_Cy3'],...
    'fontSize', fontSize,...
    'tickLength', 0.0075,...
    'figureSize', figureSize)

data_out{1} = [data.time{2}(:,end),  data.rois(i,2).timeSeries, data.rois(i,2).fit.ideal];
%% WT 
load('20230711_img02_9bp-WT-1nM_U1C-100nM_DMSO.mat')

h= figure; hold on
i = 251;
yline(data.rois(i,2).fit.components(1,2), '--', 'color', [0.5, 0.5, 0.5])
yline(data.rois(i,2).fit.components(2,2), '--', 'color', [0.5, 0.5, 0.5])
plot( data.time{end}(:,end), data.rois(i,2).timeSeries, 'color', [0.4660 0.6740 0.1880])
plot( data.time{end}(:,end), data.rois(i,2).fit.ideal, '-k')
xlim([0, 1800]);
xticks(0:300:1800); 
xlabel('Time (s)');
ylim([0, 700])
yticks(0:200:600)
% ylabel('Fluorescence (AU)')
publishFigure(h,...
    'figureName', [figurePath, 'trace_9bp_1nM_Cy3'],...
    'fontSize', fontSize,...
    'tickLength', 0.0075,...
    'figureSize', figureSize)

data_out{2} = [data.time{2}(:,end),  data.rois(i,2).timeSeries, data.rois(i,2).fit.ideal];
