close all 
clear

%% set path
figurePath = '/Users/work/Documents/Research/White_2023/figures-raw/figure01/';


%% some general stuff
fontSize = 7;
names = {'Unbound', 'Bound'};

%% Load 9bp-1A binding data
load('9bp-1A-3cy3/analysis_9bp-1A.mat')
analysis1 = analysis;

%% merge 1 nM 9bp-1A DMSO with Branaplam set
analysis = analysis1(vertcat(analysis1.U1C_nM) == 100);
nFiles  = numel(analysis);

%% simplify data
idx1 = vertcat(analysis.U1C_nM) == 100;
bootstrap = 0;
tauGuess{1} = [];
tauGuess{2} = [0.1, 100, 1000];
useFitLimits = 1;
a1 = mergeDataByField(analysis(idx1==1), 'rna_nM', 'tauGuess', tauGuess, 'fitLimits', useFitLimits, 'bootstrap', bootstrap, 'dropLastUnbound', 0);
N = numel(a1);


%% Overlay bound dwell times (CDF or PDF)
% Updated to export 
rna_nM =  vertcat(a1.field_value);
n = length(rna_nM);
raw_dwells = {};
plot_time = {};
cummulative_probability = {}; 
for i=1:2
    % unbound
    if i == 1
        [cs]= colorGradient([0.7,0.7,0.7], [0.8500 0.3250 0.0980], n);
        % bound
    else
        [cs]= colorGradient([0.7,0.7,0.7], [0 0.4470 0.7410], n);
    end
    
    h = figure;
    hold on
    ls = {};
    disp(i)
    for j = 1:length(rna_nM)
        ii = find(vertcat(a1.field_value) == rna_nM(j));
        x = a1(ii).dwells{i};
        [x_unique,~,y_cdf] = countCDF(x);
        ls{j} = [num2str(a1(ii).field_value)];
        % scatter(x_unique, y_cdf, 2, 'Marker', 'o', 'MarkerFaceColor', cs(j,:), 'MarkerEdgeColor',  cs(j,:));
        plot(x_unique, y_cdf, '-', 'color', cs(j,:));
        disp(num2str(length(x)))
        
        % store for output to excel
        raw_dwells{i,j} = x;
        plot_time{i,j} = x_unique;
        cummulative_probability{i,j} = y_cdf';
        
    end
    ylabel('Cumulative Probability')
    legend(ls,'Location', 'northwest')
    if i == 1
        xlabel('Unbound Time (s)')
        figureName = [figurePath, 'overlay_unbound_CDF'];
    elseif i==2
        xlabel('Bound Time (s)')
        figureName = [figurePath, 'overlay_bound_CDF'];
    else
        xlabel('Time To First Binding (s)')
        figureName = [figurePath, 'timeToFirstCDF'];
    end
    ylim([0,1])
    set(gca,'xscale', 'log')
    xlim([1 1e4])
    xticks([1, 1e1, 1e2, 1e3, 1e4])
    xticklabels({'10^0','10^1','10^2','10^3', '10^4'})
    publishFigure(h,...
        'figureName', figureName,...
        'fontSize', fontSize,...
        'legendFontSize', fontSize-1,...
        'tickLength', 0.03,...
        'legendTokenSize', [10,6]*0.75,...
        'figureSize', [1.75, 1.5]);
end

%% CDF overlay at a 1 nM

K = 2;
[h, plot_data ]= plotDwellsCDF(a1(K).dwells{2}, a1(K).dwellsmle{2}, 'markerSize', 3);
ylim([0,1])
set(gca,'xscale', 'log')
xlim([1 1e4])
xticks([1, 1e1, 1e2, 1e3, 1e4])
xticklabels({'10^0','10^1','10^2','10^3', '10^4'})
ylabel('Cumulative Probability')
xlabel('Bound Time (s)')
legend({['N = ' num2str(numel(a1(K).dwells{2}))], 'monoexp', 'biexp'}, 'location', 'northwest')
publishFigure(h, ...
    'figureName', [figurePath, '1nM_bound_CDF_MLE'], ...
    'legendFontSize', 6, ...
    'fontSize', 7, ...
    'tickLength', 0.03,...
    'figureSize', [1.75, 1.5]);


