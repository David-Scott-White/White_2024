%% make Figure 02 CoSMoS Example Traces
% David S. White 
% 2024-01-28

close all 
clear

%% set path
figurePath = '/Users/work/Documents/Research/White_2023/figures-raw/figure03/'; 
fontSize = 7;
figureSize = [4,0.85]

plot_data = []; 

%%  no U1C_DMSO
load('20230609_img04_9bp-1A-1nM_U1C-0nM_DMSO.mat')
roi_dmso = data.rois(398, 2);

h = figure;
hold on
yline(roi_dmso.fit.components(1,2), '--', 'color', [0.5, 0.5, 0.5])
yline(roi_dmso.fit.components(2,2), '--', 'color', [0.5, 0.5, 0.5])
plot(data.time{end}(:,end),roi_dmso.timeSeries, 'color', [0.4660 0.6740 0.1880]);
plot(data.time{end}(:,end), roi_dmso.fit.ideal, '-k')
xlabel('Time (s)')
xlim([0,1800])
xticks(0:600:1800)
% ylabel('Flourescence (au)')
yticks(0:200:800)
ylim([0,700])

plot_data = [plot_data, data.time{end}(:,end), roi_dmso.timeSeries, roi_dmso.fit.ideal]; 

publishFigure(h,...
    'figureName', [figurePath, 'exampleTrace-DMSO_noU1C'],...
    'fontSize', fontSize,...
    'tickLength', 0.02,...
    'figureSize', [2.9, 1]);

%% +Bran -U1C
load('20230612_img05_9bp-1A-1nM_U1C-0nM_Branaplam-10uM.mat')
roi_bran = data.rois(252, 2);
h = figure;
hold on
yline(roi_bran.fit.components(1,2), '--', 'color', [0.5, 0.5, 0.5])
yline(roi_bran.fit.components(2,2), '--', 'color', [0.5, 0.5, 0.5])
plot(data.time{end}(:,end),roi_bran.timeSeries, 'color', [0.4660 0.6740 0.1880]);
plot(data.time{end}(:,end), roi_bran.fit.ideal, '-k')
xlabel('Time (s)')
xlim([0,1800])
xticks(0:600:1800)
% ylabel('Flourescence (au)')
yticks(0:200:800)
ylim([0,700])

plot_data = [plot_data, data.time{end}(:,end), roi_dmso.timeSeries, roi_dmso.fit.ideal]; 

publishFigure(h,...
    'figureName', [figurePath, 'exampleTrace-Branaplam_noU1C'],...
    'fontSize', fontSize,...
    'tickLength', 0.02,...
    'figureSize', [2.9, 1]);


%% CDF overlay of 4 conditions 
