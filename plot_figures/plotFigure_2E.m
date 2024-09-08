%% make Figure 02 CoSMoS Example Traces
% David S. White 
% 2024-01-28

close all 
clear

%% set path
figurePath = '/Users/work/Documents/Research/White_2023/figures-raw/figure02/'; 
fontSize = 7;
figureSize = [4,0.85]

plot_data = []; 

%% Example Trace at 100 nM Branaplam
load('20230721_img04_9bp-1A-1nM_U1C-100nM_Branaplam-100nM')

rois = data.rois(:,2);
rois = rois(vertcat(rois.status)==1);
time_s = data.time{end}(:,end);

i = 273 % hardcoded example 215 253
h= figure; hold on
yline(rois(i).fit.components(1,2), '--', 'color', [0.5, 0.5, 0.5])
yline(rois(i).fit.components(2,2), '--', 'color', [0.5, 0.5, 0.5])
plot(time_s, rois(i).timeSeries, 'color', [0.4660 0.6740 0.1880])
plot(time_s, rois(i).fit.ideal, '-k')
xlim([0, 3600]);
xticks(0:600:3600);
xlabel('Time (s)');
ylim([0, 700])
yticks(0:200:800)
% ylabel('Fluorescence (au)')
publishFigure(h,...
    'figureName', [figurePath, '9bp-1A_1nM-100nM-Branaplam-U1C'],...
    'fontSize', 7,...
    'tickLength', 0.0075,...
    'figureSize', figureSize)

plot_data = [plot_data, time_s, rois(i).timeSeries, rois(i).fit.ideal];

%% Example Trace 10 uM Branaplam
load('20230725_img02_9bp-1A-1nM_U1C-100nM_Branaplam-10uM.mat')

rois = data.rois(:,2);
rois = rois(vertcat(rois.status)==1);
time_s = data.time{end}(:,end);

i = 201; % hardcoded example 179
h= figure; hold on
yline(rois(i).fit.components(1,2), '--', 'color', [0.5, 0.5, 0.5])
yline(rois(i).fit.components(2,2), '--', 'color', [0.5, 0.5, 0.5])
plot(time_s, rois(i).timeSeries, 'color', [0.4660 0.6740 0.1880])
plot(time_s, rois(i).fit.ideal, '-k')
xlim([0, 5400]);
xticks(0:600:5400);
xlabel('Time (s)');
ylim([0, 700])
yticks(0:200:800)
% ylabel('Fluorescence (au)')
publishFigure(h,...
    'figureName', [figurePath, '9bp-1A_1nM-10uM-Branaplam-U1C'],...
    'fontSize', 7,...
    'tickLength', 0.0075,...
    'figureSize', figureSize)

plot_data = [plot_data, time_s, rois(i).timeSeries, rois(i).fit.ideal]

