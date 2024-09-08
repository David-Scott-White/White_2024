%% Kinetic Analysis of 9bp 
% 2024-01-21
% David S. White

% requires smtoolbox.
% make sure data is in the same directory

%% Load data
close all 
clear
load('analysis_9bp.mat')

%% Set figure path
figurePath = '/Users/work/Documents/Research/White_2023/figures-raw/figureS4_9bp/';
fontSize  = 7;
figureSize = [1.75, 1.5];

%% Grab data that has saturating U1C
analysis1 = analysis(vertcat(analysis.U1C_nM)==100);

%% sum molecules 
x = vertcat(analysis1.rna_nM);
xu = unique(x); 
for i = 1:length(xu)
    xi = analysis1(vertcat(analysis1.rna_nM) == xu(i));
    disp([num2str(xu(i)), ' ', num2str(sum(vertcat(xi.nrois)))])
end

%% Kon from dwells 
% All unbound dwells
x = vertcat(analysis1.rna_nM);
y = x*0;
for i =1:length(x)
    %dtf = fitDwells(analysis1(i).dwells{1}, 1, 0, [3, 1800]);
    dtf = fitDwells(analysis1(i).timeToFirst, 1, 0, [3, 1800]);
    y(i) = dtf.monoExpTau;
end
y = 1./y(:,1);

[mdl, RMSE, R2] = linearRegression(x/1e9, y, 'intercept', 0);
kon_all = mdl.m;
ci = confint(mdl);
%kon_all
%slope_se = abs(mdl.m-ci(1,1))/1.96
h = plotLinearRegression(mdl, x/1e9, y,'xmin', 0, 'markerSize', 4, 'color', [0 0.4470 0.7410],'markerFaceColor', 'w', 'marker', '^');

%% Fit Keq
x = vertcat(analysis1.rna_nM);

close all 

% Get fits
y = zeros(length(x),2); 
yplot = cell(length(x),1); 
plotFigure = 0; 
% manual guesses
t0 = [0.001, 0.003, 0.003, 0.005, 0.003, 0.005, 0.005, 0.005];
for i = 1:length(x)
    [mdl, yplot{i}] = fitExpCurve(analysis1(i).time_s, analysis1(i).fractionbound_norm, t0(i), plotFigure);
    y(i,1) = mdl.b;
    ci = confint(mdl);
    y(i,2) = abs(mdl.b-ci(1,2));
    % title(i); 
end

%% Plot keq, all parameters free
[mdl, RMSE, R2] = linearRegression(x/1e9, y(:,1));
ci = confint(mdl); 
koff_free = mdl.b
slope_se = abs(mdl.m-ci(1,1))/1.96
intercept_se = abs(mdl.b-ci(1,2))/1.96

h = plotLinearRegression(mdl, x/1e9, y(:,1),'xmin', 0, 'markerSize', 4, 'color', [0 0.4470 0.7410],'markerFaceColor', 'w', 'marker', '^');
xlabel('[9bp] (M)')
ylabel('{\it k_{eq}} (s^{-1})')

publishFigure(h,...
    'figureName', [figurePath, 'keq_linreg'],...
    'legendFontSize', fontSize-2, ...
    'fontSize', fontSize, ...
    'tickLength', 0.03,...
    'figureSize', figureSize);

%% Linear regression of keq with a known slope
[mdl, RMSE, R2] = linearRegression(x/1e9, y(:,1), 'slope', kon_all);
koff = mdl.b
ci = confint(mdl); 
intercept_se = abs(mdl.b-ci(1))/1.96

[h, plot_data] = plotLinearRegression(mdl, x/1e9, y(:,1),'xmin', 0, 'markerSize', 4, 'color', 'r','markerFaceColor', 'w', 'marker', 'o');
xlabel('[9bp] (M)')
% ylabel('{\it k_{eq}} (s^{-1})')
ylabel('keq (s-1)')
legend('95% CI', ['R^2 = ', sprintf('%.2f', R2)], 'Data', 'location', 'northwest')
xticks([0:0.5:2]/1e9)
ax = gca;
ax.YAxis.Exponent = -2;
yticklabels({'0.2', '0', '0.2','0.4','0.6','0.8', '1'})
publishFigure(h,...
    'figureName', [figurePath, 'keq_linreg_known_kon'],...
    'legendFontSize', fontSize-2, ...
    'fontSize', fontSize, ...
    'tickLength', 0.03,...
    'figureSize', figureSize);

%% Make nice figures of equilibration traces 
color_scheme = colorGradient([0.7,0.7,0.7], [0 0.4470 0.7410], length(unique(x))+1);
% idx = [2, 1, 7, 8];
idx = [2, 4, 1, 3, 6, 8, 5, 7]
h = figure;
hold on
N = {};
colors = [[0 0.4470 0.7410]; [0.8500 0.3250 0.0980]; [0.4660 0.6740 0.1880]; [0.4940 0.1840 0.5560]];
% make replicates dashed
p = 1
for i = 1:length(idx)
    ii = idx(i); 
    N{i} = num2str(analysis1(ii).nrois);
    if rem(i, 2)
        plot(analysis1(ii).time_s, analysis1(ii).fractionbound_norm, '-', 'color', colors(p,:))
    else
        plot(analysis1(ii).time_s, analysis1(ii).fractionbound_norm, '--', 'color', colors(p,:))
        p = p+1
    end
end
for i = 1:length(idx)
    ii = idx(i); 
    plot(analysis1(ii).time_s, yplot{ii}, '-k')
end

xlabel('Time (s)')
ylabel('Normalized Fraction Bound')
legend({['0.25 nM, N = ', N{1}],...
    ['0.25 nM, N = ', N{2}],...
    ['0.5 nM, N = ', N{3}],...
    ['0.5 nM, N = ', N{4}],...
    ['1 nM, N = ', N{5}],...
    ['1 nM, N = ', N{6}],...
    ['2 nM, N = ', N{7}],...
    ['2 nM, N = ', N{8}]}, 'location', 'southeast')
xlim([0, 1800])
xticks(0:600:1800)
publishFigure(h,...
    'figureName', [figurePath, 'fractionbound_keq'],...
    'legendFontSize', fontSize, ...
    'fontSize', fontSize, ...
    'tickLength', 0.03,...
    'figureSize', [3,2.75]);

%% Bound dwell times, merge concentrations 
% rna_nM = unique(x); 
% for i = 1:length(rna_nM)
%     a = analysis1(vertcat(analysis1.rna_nM) == rna_nM(i)); 
%     dwells = []; 
%     for j = 1:length(a)
%         dwells = [dwells; a(j).dwells{2}];
%     end
%     dtf = fitDwells(dwells, 1, 0, []); 
%     
%     h = plotDwells(dwells, dtf, 'legendLocation', 'northwest'); 
%     delete(h.Children(end).Children(1))
%     legend({['N = ', num2str(dtf.nDwells)]})
%     xlim([0, 1800])
%     xticks(0:600:1800)
%     xlabel('Bound Time (s)')
%     publishFigure(h,...
%     'figureName', [figurePath, ['dwells_', num2str(rna_nM(i))]],...
%     'legendFontSize', fontSize-1, ...
%     'fontSize', fontSize, ...
%     'tickLength', 0.03,...
%     'figureSize', figureSize);
% end
% 
% %% Unbound dwell times, merge concentrations 
% rna_nM = unique(x); 
% for i = 1:length(rna_nM)
%     a = analysis1(vertcat(analysis1.rna_nM) == rna_nM(i)); 
%     dwells = []; 
%     for j = 1:length(a)
%         dwells = [dwells; a(j).dwells{1}];
%     end
%     dtf = fitDwells(dwells, 1, 0); 
%     
%     h = plotDwells(dwells, dtf, 0, [3,1800], 'legendLocation', 'northeast', 'faceColor', [1,1,1]); 
%     legend({['N = ', num2str(dtf.nDwells)], 'monoexp'})
%     xlim([0, 1800])
%     xticks(0:600:1800)
%     xlabel('Unbound Time (s)')
%     publishFigure(h,...
%     'figureName', [figurePath, ['dwells_unbound_', num2str(rna_nM(i))]],...
%     'legendFontSize', fontSize-1, ...
%     'fontSize', fontSize, ...
%     'tickLength', 0.03,...
%     'figureSize', figureSize);
% end

%% Time to first event dwell times 
% rna_nM = unique(x); 
% for i = 1:length(rna_nM)
%     a = analysis1(vertcat(analysis1.rna_nM) == rna_nM(i)); 
%     dwells = []; 
%     for j = 1:length(a)
%         dwells = [dwells; a(j).timeToFirst];
%     end
%     dtf = fitDwells(dwells, 1, 0); 
%     
%     h = plotDwells(dwells, dtf, 0, [3, 1800], 'legendLocation', 'northeast', 'faceColor', [1,1,1]); 
%     xlim([0, 1800])
%     xticks(0:300:1800)
%     xlabel('Time To First Binding (s)')
%     publishFigure(h,...
%     'figureName', [figurePath, ['dwells_unbound_', num2str(rna_nM(i))]],...
%     'legendFontSize', fontSize-1, ...
%     'fontSize', fontSize, ...
%     'tickLength', 0.03,...
%     'figureSize', figureSize);
% end

%% plot raster
h = plotRaster(analysis(6).events, 'frameRate_s', 3);
xticks(0:600:1800)
publishFigure(h,...
    'figureName', [figurePath, 'raster_1nM'],...
    'legendFontSize', fontSize, ...
    'fontSize', fontSize, ...
    'tickLength', 0.03,...
    'figureSize', [3,2.75]);

