%% Make Figure 04 
% David S. White 
% 2023-09-19

% still figuring this one out. the goal should be to show something about
% how the ka and koff are conserved across non -1A bulged oligos 

%% Clean
close all 
clear

%% put stuff here
figurePath = '/Users/work/Documents/Research/White_2023/figures-raw/figure05/';

%% Some data. lets load it. Make sure single-molecule data is in filepath... too lazy to write it all in for loading otherwise
% no analysis file for mimic since there was no binding. load otherwise if
% needed
%names = {'+2C', '-1C', '6-Exon', '6-Intron', '+1C'};
%files = {'analysis_9bp2C.mat', 'analysis_9bp-1C.mat', 'analysis_6bp-Exon.mat', 'analysis_6bp-Intron.mat', 'analysis_9bp1C.mat'};

names = {'-1C', '+1C', '+2C', '6bp-exon', '6-Intron'};
files = {'analysis_9bp-1C.mat', 'analysis_9bp1C.mat', 'analysis_9bp2C.mat', 'analysis_6bp-Exon.mat', 'analysis_6bp-Intron.mat'}

nFiles = length(files); 
A = cell(nFiles,1); 
for i = 1:nFiles
    load([files{i}]); 
    A{i} = analysis; 
end

%% now what to plot... 
% Maybe the intersting things it the WT and +2C, then the handedness by
% 6-exon & 6-intron. Move -1C to supplement. Move +1C to supplement with
% the mimic 

%% Compare data. Option 1, plot CDF of unbound and bound times at +/- U1C
A0 = cell(5,1); 
for i = 1:length(A)
    idx1 = find(vertcat(A{i}.U1C_nM) == 0);
    idx2 = find(vertcat(A{i}.U1C_nM) == 100);
    idx3 = find(vertcat(A{i}.cpd_uM) == 0);
    idx1 = intersect(idx1, idx3); 
    idx2 = intersect(idx2, idx3);
    a0 = mergeDataByField(A{i}(idx1), 'rna_nM', 'tauGuess',[], 'fitLimits', 1, 'bootstrap', 0);
    a1 = mergeDataByField(A{i}(idx2), 'rna_nM', 'tauGuess',[], 'fitLimits', 1, 'bootstrap', 0);
    A0{i} = [a0; a1];
end

%% now lets make them as violin plots
plot_data = {}; 
for i = 1:length(A0)
    display(names{i})
    ii = 1;
    for k = 1:2
        X = [];
        Y = [];
        str = cell(2,1); 
        for j = 1:2
            x = A0{i}(j).dwells{k};
            y = x*0 + (j-1)*100;
            X = [X;x];
            Y = [Y;y];
            str{j} = [sprintf('%s%i\n%0.1f', 'N = ', numel(x), A0{i}(j).dwellsmle{k}.monoExpTau), ' s'];
            plot_data{i, ii} = x;
            ii = ii +1;
        end
       
        h=figure;
        vp = violinplot(log10(X), Y, 'ShowData', false, 'EdgeColor', [0,0,0],...
            'MedianColor', [1,1,1], 'BoxColor',  [0,0,0], 'ViolinAlpha', 0.75);
        ylim([0,4])
        yticks(0:1:4)
        yticklabels({'10^{0}','10^{1}','10^{2}','10^{3}', '10^{4}'})
        if k == 1
            ylabel('Unbound Time (s)')
            figname = 'unbound_violin_';
        else
            ylabel('Bound Time (s)')
            figname = 'bound_violin_';
        end
        text(1,4, str{1}, 'FontName', 'Arial','HorizontalAlignment', 'center', 'FontSize', 6);
        text(2,4, str{2}, 'FontName', 'Arial','HorizontalAlignment', 'center', 'FontSize', 6);
        text(1.5, 4.7, names{i},  'FontName', 'Arial','HorizontalAlignment', 'center', 'FontSize', 6, 'FontWeight', 'bold');
        xlabel('[U1C] (nM)')
        xlim([0.5, 2.5])
        publishFigure(h, ...
            'figureName', [figurePath, figname, names{i}], ...
            'legendFontSize', 6, ...
            'fontSize', 7, ...
            'tickLength', 0.03,...
            'figureSize', [1.5, 1.5]);
    end
    
end



