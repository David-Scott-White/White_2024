%% make Figure 02 CoSMoS Branaplam Analysis
% David S. White 
% 2024-01-28

close all 
clear

%% set path
figurePath = '/Users/work/Documents/Research/White_2023/figures-raw/figure03/';


%% some general stuff
fontSize = 7;
names = {'Unbound', 'Bound'};

%% Load 9bp-1A binding data
load('9bp-1A-3cy3/analysis_9bp-1A.mat')
analysis1 = analysis;

%% Branaplam binding set
load('analysis_9bp-1A_Branaplam.mat')
analysis2 = analysis;

%% merge 1 nM 9bp-1A DMSO with Branaplam set
idx1 = find(vertcat(analysis1.rna_nM) == 1);
idx2 = find(vertcat(analysis1.U1C_nM) == 100);
idx3 = intersect(idx1, idx2);
analysis = [analysis1(idx3); analysis2];
nFiles  = numel(analysis);

%% simplify data
idx1 = vertcat(analysis.U1C_nM) == 100;
bootstrap = 0;
tauGuess{1} = [];
tauGuess{2} = [0.1, 100, 1000];
useFitLimits = 1;
a1 = mergeDataByField(analysis(idx1==1), 'cpd_uM', 'tauGuess', tauGuess, 'fitLimits', useFitLimits, 'bootstrap', bootstrap, 'dropLastUnbound', 0);


%% Load in -U1C 1- uM data
idx1 = find(vertcat(analysis.cpd_uM) == 10);
bootstrap = 0;
tauGuess{1} = [];
tauGuess{2} = [0.1, 100, 1000];
useFitLimits = 1;
a2 = mergeDataByField(analysis(idx1), 'U1C_nM', 'tauGuess', tauGuess, 'fitLimits', useFitLimits, 'bootstrap', bootstrap);
N = numel(a2);

%% merge data into 4 categories 
% its messy but hey it works. 
A00 = [analysis1; analysis2]; 
idx1 = find(vertcat(A00.cpd_uM) == 0);
idx2 = find(vertcat(A00.rna_nM) == 1);
idx3 = find(vertcat(A00.cpd_uM) == 10);

B1 = A00(intersect(idx1, idx2));
b1 = find(vertcat(B1.U1C_nM) == 0);
b2 = find(vertcat(B1.U1C_nM) == 100);

a1 = mergeDataByField(B1(b1), 'U1C_nM', 'tauGuess', {[], [0.5, 50, 100]}, 'fitLimits', 1); 
a2 = mergeDataByField(B1(b2), 'U1C_nM', 'tauGuess', {[], [0.5, 50, 300]}, 'fitLimits', 1); 

B2 = A00(intersect(idx2, idx3));
b3 = find(vertcat(B2.U1C_nM) == 0);
b4 = find(vertcat(B2.U1C_nM) == 100);

a3 = mergeDataByField(B2(b3), 'U1C_nM', 'tauGuess', {[], [0.5, 50, 50]}, 'fitLimits', 1); 
a4 = mergeDataByField(B2(b4), 'U1C_nM', 'tauGuess', {[], [0.5, 50, 1000]}, 'fitLimits', 1); 

% merge 
A00 = [];
A00 = [a1; a3; a2; a4];

%% Plot violin plot of all four options
str = cell(4,1);
X = [];
Y = [];
plot_data = {}; 

for i = 1:4
    x = A00(i).dwells{2};
    
    plot_data{i} = x; 
    
    y = x*0 + i;
    X = [X;x];
    Y = [Y;y];
    %str{i} = [sprintf('%s%i\n%0.1f', 'N = ', numel(x), A00(i).dwellsmle{2}.monoExpTau), ' s'];
    str{i} = [num2str(numel(x))];
end

h=figure;
vp = violinplot(log10(X), Y, 'ShowData', false, 'EdgeColor', [0,0,0],...
    'MedianColor', [1,1,1], 'BoxColor',  [0,0,0], 'ViolinAlpha', 0.75);
hold on
% yline(log10(mean( A00(i).dwells{2})),'-r')
ylim([0,4])
yticks(0:1:5)
yticklabels({'10^{0}','10^{1}','10^{2}','10^{3}', '10^{4}'})

ylabel('Bound Time (s)')
figname = 'bound_violin_';

for i = 1:4
    text(i,4, str{i}, 'FontName', 'Arial','HorizontalAlignment', 'center', 'FontSize', 5);
end
% xlabel('[U1C] (nM)')
xlim([0.5, 4.5])
publishFigure(h, ...
    'figureName', [figurePath, figname], ...
    'legendFontSize', 6, ...
    'fontSize', 7, ...
    'tickLength', 0.03,...
    'figureSize', [1.75, 1.5]);

