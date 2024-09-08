%% Plot Figure 1. SPR data
% David S. White 
% 2023-12-29 

% add White_2023 to file path

%% set figure path
figurePath = '/Users/work/Documents/Research/White_2023/figures-raw/figure01/';

figureSize = [1.75, 1.25];
fontSize = 7;


%% Load DMSO data(.csv from Bryan). 
% concentrations are hard coded from reading the sheet. Data in duplicate
rna_DMSO_nM = [0.01953125, 0.390625, 0.78125, 1.5625, 3.125, 6.25, 12.5, 25, 50, 100];

spr_1A_DMSO = readmatrix('Fig 1 and S2 - 1A_bulge - raw data.xlsx');
spr_11_DMSO = readmatrix('Fig 1 and S2 - match - raw data.xlsx');

spr_1A_DMSO_fit = readmatrix('Fig 1 and S2 - 1A_bulge - 1-state fits.xlsx');
spr_11_DMSO_fit = readmatrix('Fig 1 and S2 - match - 1-state fits.xlsx');

%% plot -1A
time_s = spr_1A_DMSO(:,1); 
plot_idx = 2:4:40;
cs = colorGradient([0.7,0.7,0.7], [0.4940 0.1840 0.5560], numel(plot_idx));
h=figure; hold on
plot_data_11bp1A  = time_s
for i = 1:length(plot_idx)
    plot(time_s, spr_1A_DMSO(:,plot_idx(i)), 'color', cs(i,:))
    plot_data_11bp1A = [plot_data_11bp1A, spr_1A_DMSO(:,plot_idx(i))];
end
xlabel('Time (s)')
ylabel('Response Units (RU)')
xlim([-40, 800])
ylim([0,30])
publishFigure(h,...
    'figureName', [figurePath, '11bp-1A_DMSO_spr'],...
    'legendFontSize', fontSize-2, ...
    'fontSize', fontSize, ...
    'tickLength', 0.03,...
    'figureSize', figureSize);



% Plot with single site fit
% for i = 1:length(plot_idx)
%     plot(time_s, spr_1A_DMSO_fit(:,plot_idx(i)), '-k')
% end
% xlabel('Time (s)')
% ylabel('Response Units (RU)')
% xlim([-40, 800])
% ylim([0,30])
% publishFigure(h,...
%     'figureName', [figurePath, '11bp-1A_DMSO_spr'],...
%     'legendFontSize', fontSize-2, ...
%     'fontSize', fontSize, ...
%     'tickLength', 0.03,...
%     'figureSize', figureSize);


%% plot 11bp 
cs = colorGradient([0.7,0.7,0.7], [0.4940 0.1840 0.5560], numel(plot_idx));
h=figure; hold on
plot_data_11bp = {}; 
plot_data_11bp = time_s;
for i = 1:length(plot_idx)
    plot(time_s, spr_11_DMSO(:,plot_idx(i)), 'color', cs(i,:))
    plot_data_11bp = [plot_data_11bp, spr_11_DMSO(:,plot_idx(i))];
end
xlabel('Time (s)')
ylabel('Response Units (RU)')
xlim([-40, 800])     
ylim([0,30])
publishFigure(h,...
    'figureName', [figurePath, '11bp_DMSO_spr'],...
    'legendFontSize', fontSize-2, ...
    'fontSize', fontSize, ...
    'tickLength', 0.03,...
    'figureSize', figureSize);


%% Plot with fits, -1A
C = [colororder;colororder];
time_s = spr_1A_DMSO(:,1); 
h=figure; hold on
p = 1;
for i = 2:4:40
    plot(time_s, spr_1A_DMSO(:,i), 'Color', C(p,:))
    p = p+1; 
end
p = 1;
for i = 4:4:40
    plot(time_s, spr_1A_DMSO(:,i), 'Color', C(p,:))
    p = p+1;
end
for i = 2:2:40
   plot(time_s, spr_1A_DMSO_fit(:,i), '-k')
end
xlabel('Time (s)')
ylabel('Response Units (RU)')
xlim([-40, 800])
ylim([0,30])
publishFigure(h,...
    'figureName', [figurePath, '11bp-1A_DMSO_spr_fits'],...
    'legendFontSize', fontSize-1, ...
    'fontSize', fontSize+1, ...
    'tickLength', 0.03,...
    'figureSize', figureSize*1.5);

h=figure; hold on
p = 1;
for i = 2:4:40
    plot(time_s, spr_11_DMSO(:,i), 'Color', C(p,:))
    p = p+1; 
end
p = 1;
for i = 4:4:40
    plot(time_s, spr_11_DMSO(:,i), 'Color', C(p,:))
    p = p+1;
end
for i = 2:2:40
   plot(time_s, spr_11_DMSO_fit(:,i), '-k')
end
xlabel('Time (s)')
ylabel('Response Units (RU)')
xlim([-40, 800])
ylim([0,30])
publishFigure(h,...
    'figureName', [figurePath, '11bp_DMSO_spr_fits'],...
    'legendFontSize', fontSize-1, ...
    'fontSize', fontSize+1, ...
    'tickLength', 0.03,...
    'figureSize', figureSize*1.5);


