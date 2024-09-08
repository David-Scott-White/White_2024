
%% Clean
close all 
clear

%% put stuff here
figurePath = '/Users/work/Documents/Research/White_2023/figures-raw/figure05/';

%% Some data. lets load it. Make sure single-molecule data is in filepath... too lazy to write it all in for loading otherwise
% no analysis file for mimic since there was no binding. load otherwise if
% needed
names = {'WT', '9bp-1A'};
files = {'analysis_9bp.mat', 'analysis_9bp-1A.mat'}
nFiles = length(files); 
A = cell(nFiles,1); 
for i = 1:nFiles
    load([files{i}]); 
    A{i} = analysis; 
end

%% Linear regression of kon WT
idx1 = vertcat(A{1}.U1C_nM) == 0; 
idx2 = vertcat(A{1}.U1C_nM) == 100; 

disp('> 9bp')
a0  = A{1}(idx1);
x = vertcat(a0.rna_nM);
y = x*0;
for i =1:length(x)
    dtf = fitDwells(a0(i).dwells{1}, 1, 0, [3,1800]);
    y(i) = dtf.monoExpTau;
end
y = 1./y(:,1);
disp('>> 0 nM U1C')
[mdl, RMSE, R2] = linearRegression(x, y, 'intercept', 0);
[h, plot_data_WT_0] = plotLinearRegression(mdl, x, y,'xmin', 0, 'markerSize', 4, 'color', [0.1,0.1,0.1],'markerFaceColor', 'w');

a1  = A{1}(idx2);
x = vertcat(a1.rna_nM);
y = x*0;
for i =1:length(x)
    dtf = fitDwells(a1(i).dwells{1}, 1, 0, [3,1800]);
    y(i) = dtf.monoExpTau;
end
y = 1./y(:,1);
disp('>> 100 nM U1C')
[mdl, RMSE, R2] = linearRegression(x, y, 'intercept', 0);
[~, plot_data_WT_100] = plotLinearRegression(mdl, x, y, 'parent', h,'color',[0 0.4470 0.7410], 'xmin', 0, 'marker', 'd', 'markerFaceColor', 'w', 'markerSize', 4);
xlim([0,4])
xticks(0:1:4)
xlabel('[9bp] (nM)')
ax = gca;
ax.YAxis.Exponent = -3;
ylabel('{\it k_{apparent}} (s^{-1})')
ylim([0,0.01])

publishFigure(h,...
    'figureName', [figurePath, 'WT_kon_linear-regression'],...
    'fontSize', 7,...
    'legendFontSize', 6,...
    'tickLength', 0.03,...
    'legendTokenSize', [10,6],...
    'figureSize', [1.75, 1.5]);


%% Linear regression of kon -1A
idx1 = vertcat(A{2}.U1C_nM) == 0; 
idx2 = vertcat(A{2}.U1C_nM) == 100; 

disp(' ')
disp('> 9bp-1A')
a0  = A{2}(idx1);
x = vertcat(a0.rna_nM);
y = x*0;
for i =1:length(x)
    dtf = fitDwells(a0(i).dwells{1}, 1, 0, [3,1800]);
    y(i) = dtf.monoExpTau;
end
y = 1./y(:,1);

disp(' ')
disp('>> 0 nM U1C')
[mdl, RMSE, R2] = linearRegression(x, y, 'intercept', 0);
[h, plot_data_1A_0] = plotLinearRegression(mdl, x, y,'xmin', 0, 'markerSize', 4, 'color', [0.1,0.1,0.1],'markerFaceColor', 'w');

a1  = A{2}(idx2);
x = vertcat(a1.rna_nM);
y = x*0;
for i =1:length(x)
    dtf = fitDwells(a1(i).dwells{1}, 1, 0, [3,1800]);
    y(i) = dtf.monoExpTau;
end
y = 1./y(:,1);

disp(' ')
disp('>> 100 nM U1C')
[mdl, RMSE, R2] = linearRegression(x, y, 'intercept', 0);
[~, plot_data_1A_100] = plotLinearRegression(mdl, x, y, 'parent', h,'color', [0.4660 0.6740 0.1880], 'xmin', 0, 'marker', 'd', 'markerFaceColor', 'w', 'markerSize', 4);
xlim([0,4])
xticks(0:1:4)
xlabel('[9bp-1A] (nM)')
ax = gca;
ax.YAxis.Exponent = -3;
ylabel('{\it k_{apparent}} (s^{-1})')
ylim([0,0.015])

publishFigure(h,...
    'figureName', [figurePath, '9bp1A_kon_linear-regression'],...
    'fontSize', 7,...
    'legendFontSize', 6,...
    'tickLength', 0.03,...
    'legendTokenSize', [10,6],...
    'figureSize', [1.75, 1.5]);