%% Make Figure 05 
% Plot bound dwell comparion of alterantive splicing mutants +/- U1C, +/-
% Branpalam

%% Figure path

figurePath = '/Users/work/Documents/Research/White_2023/figures-raw/figure03/';
fontSize = 7

%%  Load data, make sure single-molecule-data is in file path
names = {'SMN2', 'FOXM1', 'HTT', 'SF3B3'};
files = {'analysis_SMN2.mat', 'analysis_FOXM1.mat', 'analysis_HTT.mat', 'analysis_SF3B3.mat'};
nFiles = length(files); 
A = cell(nFiles,1); 
for i = 1:nFiles
    load(files{i})
    % index
    idx1 = find(vertcat(analysis.U1C_nM) == 0);
    idx2 = find(vertcat(analysis.U1C_nM) == 100);
    idx3 = find(vertcat(analysis.cpd_uM) == 0);
    idx4 = find(vertcat(analysis.cpd_uM) == 10);
    
    % no U1C, DMSO
    a1 = mergeDataByField(analysis(intersect(idx1, idx3)), 'cpd_uM', 'tauGuess',[], 'fitLimits', 1, 'bootstrap', 0);
    % no U1C, +Branaplam
    a2 = mergeDataByField(analysis(intersect(idx1, idx4)), 'cpd_uM', 'tauGuess',[], 'fitLimits', 1, 'bootstrap', 0);
    % +U1C, DMSO
    a3 = mergeDataByField(analysis(intersect(idx2, idx3)), 'cpd_uM', 'tauGuess',[], 'fitLimits', 1, 'bootstrap', 0);
    % +U1C, +Branaplam
    a4 = mergeDataByField(analysis(intersect(idx2, idx4)), 'cpd_uM', 'tauGuess',[], 'fitLimits', 1, 'bootstrap', 0);
    
    A{i} = [a1;a2;a3;a4];
end

%% Plot Violin plots
plot_data_D = {}; 
for i = 1:4
    X = [];
    Y = [];
    N = {}; 
    for j = 1:4
        x = A{i}(j).dwells{2};
        if length(x) < 50
            x = 0; 
        end
        plot_data_D{i,j} = x; 
        display(names{i})
        N{j} = num2str(length(x)); 
        % figure; plotDwells(A{i}(j).dwells{2}, A{i}(j).dwellsmle{2})
        y = x*0+j;
        X = [X; x];
        Y = [Y; y];
    end
    
    h=figure;
    vp = violinplot(log10(X), Y, 'ShowData', false, 'EdgeColor', [0,0,0],...
        'MedianColor', [1,1,1], 'BoxColor',  [0,0,0], 'ViolinAlpha', 0.75);
    ylim([0,3.5])
    yticks(0:1:4)
    yticklabels({'10^{0}','10^{1}','10^{2}','10^{3}', '10^{4}'})
    %title(names{i})
    ylabel('Bound Time (s)')
    xlim([0.5, 4.5])
    for j = 1:4
        text(j,3.4, N{j}, 'FontName', 'Arial','HorizontalAlignment', 'center', 'FontSize', 6);
    end

    publishFigure(h, ...
        'figureName', [figurePath, 'bound_violin_', names{i}], ...
        'legendFontSize', 8, ...
        'fontSize', 7, ...
        'tickLength', 0.03,...
        'figureSize', [2, 1.75]);
end

%% make bar graph of color gradient for AI use
% ncolor = 4;
% [grad,im]=colorGradient([1,1,1],[0.4940 0.1840 0.5560],ncolor);
% h = figure; hold on
% for i = 1:ncolor
%     bar(i, 10, 'FaceColor', grad(i,:))
% end
% publishFigure(h, ...
%     'figureName', [figurePath, '_purple'], ...
%     'legendFontSize', 8, ...
%     'fontSize', 7, ...
%     'tickLength', 0.03,...
%     'figureSize', [2, 1.75]);
% 
% 
% [grad,im]=colorGradient([1,1,1],[0 0.4470 0.7410],ncolor);
% h = figure; hold on
% for i = 1:ncolor
%     bar(i, 10, 'FaceColor', grad(i,:))
% end
% publishFigure(h, ...
%     'figureName', [figurePath, '_blue'], ...
%     'legendFontSize', 8, ...
%     'fontSize', 7, ...
%     'tickLength', 0.03,...
%     'figureSize', [2, 1.75]);
% 
% [grad,im]=colorGradient([1,1,1],[0 0.4470 0.7410],ncolor);
% h = figure; hold on
% for i = 1:ncolor
%     bar(i, 10, 'FaceColor', grad(i,:))
% end
% publishFigure(h, ...
%     'figureName', [figurePath, '_blue'], ...
%     'legendFontSize', 8, ...
%     'fontSize', 7, ...
%     'tickLength', 0.03,...
%     'figureSize', [2, 1.75]);
% 
% % RED
% [grad,im]=colorGradient([1,1,1],[0.6350 0.0780 0.1840],ncolor);
% h = figure; hold on
% for i = 1:ncolor
%     bar(i, 10, 'FaceColor', grad(i,:))
% end
% publishFigure(h, ...
%     'figureName', [figurePath, '_red'], ...
%     'legendFontSize', 8, ...
%     'fontSize', 7, ...
%     'tickLength', 0.03,...
%     'figureSize', [2, 1.75]);

