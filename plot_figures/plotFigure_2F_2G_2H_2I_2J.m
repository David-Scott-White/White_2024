%% make Figure 02 CoSMoS Branaplam Analysis
% David S. White 
% 2024-01-28

close all 
clear

%% set path
figurePath = '/Users/work/Documents/Research/White_2023/figures-raw/figure02/';

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
N = numel(a1);

%% Plot fraction bound as function of cpd_uM
A0 = analysis(idx1==1);
x = vertcat(A0.cpd_uM); 
y = vertcat(A0.fractionbound_fit);


% fit EC50 (with hill slope set to 1)
func = @(beta,x) beta(1) + ((beta(2)-beta(1)) ./ (1+beta(3)./x));
minResponse=min(y(:,1));
maxResponse=max(y(:,1));
EC50_guess = 3;
[coeffs,r,J]=nlinfit(x(:),y(:,1),func,[minResponse, maxResponse EC50_guess]);
CI = nlparci(coeffs,r,'jacobian',J);
se = (coeffs(end) - CI(3,1)) / 1.96;

disp(['EC50: ', num2str(coeffs(3)), ' ± ', num2str(se)])

% EC50 Log break
xrange = logspace(-6, 1, 500);
xrange = xrange(:);
[Ypred,delta] = nlpredci(func,xrange,coeffs,r,J);
Ypred_ci = delta*[-1 1]+Ypred;
fig3 = figure;
hold on
xrange = log10(xrange);
patch([xrange' fliplr(xrange')], [Ypred_ci(:,2)' fliplr(Ypred_ci(:,1)')], '',...
    'FaceColor',[1 0 0], 'EdgeColor','none', 'FaceAlpha', 0.2);
plot(xrange, Ypred, '-', 'Color', [1 0 0])

% mean fraction bound, error = standard error
% x(1) = 1e-6; % doesnt look like the prediction changes below this value as it reaches zero
x(x==0) = 1e-6;
% errorbar(log10(x),y (:,1), y (:,3), y (:,3), 'ok', 'MarkerSize', 5, 'CapSize', 5, 'MarkerFaceColor', 'w')
scatter(log10(x),y (:,1), 15, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w')
ylabel('Fraction Bound')
xlabel('[Branaplam] (µM)')
ylim([0.4, 1])

set(gca,'tickdir', 'out', 'fontSize', 7,'TickLength', [0.03, 0.03],'FontName', 'Arial')
fig3.Units = 'inches';
fig3.Position = [0, 0, 1.75, 1.5];
xlim([-6.25, 1])
xticks(-6:1:-1)
xticklabels({'0.1', '1', '10'})
set(gcf,'color','w');
set(gca, 'Layer', 'top', 'XColor', 'k', 'YColor', 'k', 'TickLength', [0.03, 0.03], 'FontName', 'Arial')
h = breakxaxis([-5.75, -1.25]);
h.leftAxes.XTick = -6;
h.leftAxes.XTickLabel = '0';
h.rightAxes.XTick = [-1,0,1];
exportgraphics(gcf,[figurePath, '9bp-1A_Branaplam-EC50_log_break.pdf'], 'ContentType','vector');

plot_data_F_raw = [x,y(:,1)];
plot_data_F_fit = [xrange, Ypred, Ypred_ci];


%% Overlay bound dwell times (CDF or PDF)
branaplam_uM =  vertcat(a1.field_value);
branaplam_uM = branaplam_uM([1:1:6]);
n = length(branaplam_uM);
plot_data_G = {};
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
    for j = 1:length(branaplam_uM)
        ii = find(vertcat(a1.field_value) == branaplam_uM(j));
        x = a1(ii).dwells{i};
        [x_unique,~,y_cdf] = countCDF(x);
        % ls{j} = [num2str(a1(ii).field_value), ' µM'];
        ls{j} = [num2str(a1(ii).field_value)];
        scatter(x_unique, y_cdf, 2, 'Marker', 'o', 'MarkerFaceColor', cs(j,:), 'MarkerEdgeColor',  cs(j,:));
        disp(num2str(length(x)))
        
        plot_data_G{i,j} = [x_unique(:), y_cdf(:)];
        
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


%% Bound mle tau 
x_uM = vertcat(a1.field_value);
x_uM(1)=1e-6;
cs = colorGradient([1,1,1], [0 0.4470 0.7410], 3);
boundmle = vertcat(a1.boundmle);
boundmle_se = vertcat(a1.boundmle_se);

h1 = figure; hold on
plot(log10(x_uM), boundmle(:,3), '--', 'Color', [0.3010 0.7450 0.9330])
scatter(log10(x_uM), boundmle(:,3), 15, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.3010 0.7450 0.9330])
errorbar(log10(x_uM), boundmle(:,3), boundmle_se(:,3), boundmle_se(:,3), 'LineStyle', 'None', 'Color', 'k', 'MarkerFaceColor', [0.3010 0.7450 0.9330], 'MarkerSize', 0.1, 'CapSize', 0, 'Marker', 'o');


plot(log10(x_uM), boundmle(:,4), '--', 'Color', [0.8500 0.3250 0.0980])
scatter(log10(x_uM), boundmle(:,4), 15, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.8500 0.3250 0.0980])
errorbar(log10(x_uM), boundmle(:,4), boundmle_se(:,4), boundmle_se(:,4), 'Color', 'k', 'LineStyle', 'None', 'MarkerFaceColor', [0.8500 0.3250 0.0980], 'MarkerSize', 0.1, 'CapSize', 0, 'Marker', 'o');
legend({'{\it tau_{B}^1}','{\it tau_{B}^2}'}, 'location','best')
print(gcf, [figurePath, '9bp-1A_U1C_bound-tau_MLE'], '-dpdf')


% previous label but it messes with axes size.. fix in Adobe
% ylabel('{\it k_{off}} (s^{-1})')
ylabel('tau (s)')
xlabel('[Branaplam] (µM)')
%ylim([1e-4, 1e-1])
set(gca,'yscale','log')
ax = gca;
%ax.YAxis.Exponent = -2;
%yticklabels({'0.4', '0.3', '0.2', '0.1'})

set(gca,'tickdir', 'out', 'FontSize', 7, 'FontName', 'Arial',...
    'box', 'off', 'Layer', 'top', 'XColor', 'k', 'YColor', 'k','TickLength', [0.03, 0.03])
h1.Units = 'inches';
h1.Position = [0, 0, 1.75, 1.5];
xlabel('[Branaplam] (µM)')
% xlim([-6.25, 2.25])
xlim([-6.25, 1])
xticks(-6:1:-1)
xticklabels({'10^{-1}', '10^{0}', '10^{1}', '10^{2}'})
set(gcf,'color','w');

legend({'{\it k_{off}^1}','{\it k_{off}^2}'}, 'location','best')
print(gcf, [figurePath, '9bp-1A_branaplam_bound-tau_MLE_legend'], '-dpdf')

h = breakxaxis([-5.75, -1.25]);
h.leftAxes.XTick = -6;
h.leftAxes.XTickLabel = '0';
h.rightAxes.XTick = [-1,0,1,2];
legend({'{\it k_{off}^1}','{\it k_{off}^2}'}, 'location','best')
legend('off')
print(gcf, [figurePath, '9bp-1A_branaplam_bound-tau_MLE'], '-dpdf')


plot_data_I = [x_uM, boundmle(:,3),  boundmle_se(:,3), boundmle(:,4), boundmle_se(:,4)];

%% Bound MLE Amplitude
h1 = figure; hold on
boundmle = vertcat(a1.boundmle);
boundmle_se = vertcat(a1.boundmle_se);

plot(log10(x_uM),boundmle(:,2), '-', 'Color', [0.3010 0.7450 0.9330])
scatter(log10(x_uM),boundmle(:,2), 15, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.3010 0.7450 0.9330])
errorbar(log10(x_uM),boundmle(:,2), boundmle_se(:,2), boundmle_se(:,2), 'LineStyle', 'None',  'Color', 'k', 'MarkerFaceColor', [0.3010 0.7450 0.9330], 'MarkerSize', 0.1, 'CapSize', 0, 'Marker', 'o');

plot(log10(x_uM),1-boundmle(:,2), '-', 'Color', [0.8500 0.3250 0.0980])
scatter(log10(x_uM), 1-boundmle(:,2), 15, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.8500 0.3250 0.0980])
errorbar(log10(x_uM), 1-boundmle(:,2), boundmle_se(:,2),boundmle_se(:,2), 'LineStyle', 'None', 'Color', 'k', 'MarkerFaceColor', [0.8500 0.3250 0.0980], 'MarkerSize', 0.1, 'CapSize', 0, 'Marker', 'o');

legend({'{\it k_{off}^1}','{\it k_{off}^2}'}, 'location','best');
print(gcf, [figurePath, '9bp-1A_U1C_bound-koff_MLE'], '-dpdf');

ylabel('Amplitude')
xlabel('[Branaplam] (µM)')
ylim([0, 1.1])
%legend({'{\it k_{off}^1}','{\it k_{off}^2}'}, 'location','best')
set(gca,'tickdir', 'out', 'FontSize', 7, 'FontName', 'Arial',...
    'box', 'off', 'Layer', 'top', 'XColor', 'k', 'YColor', 'k','TickLength', [0.03, 0.03]);
h1.Units = 'inches';
h1.Position = [0, 0, 1.75, 1.5];
xlabel('[Branaplam] (µM)')
xlim([-6.25, 1]);
xticks(-6:1:-1);
xticklabels({'10^{-1}', '10^{0}', '10^{1}', '10^{2}'})
set(gcf,'color','w');
h = breakxaxis([-5.75, -1.25]);
h.leftAxes.XTick = -6;
h.leftAxes.XTickLabel = '0';
h.rightAxes.XTick = [-1,0,1,2];
print(gcf, [figurePath, '9bp-1A_branaplam_bound-Amp_MLE'], '-dpdf')

plot_data_I = [plot_data_I, 1-boundmle(:,2),  boundmle_se(:,2)];
    
%% Bound dwell fit with hill 
x_uM = vertcat(a1.field_value);
x_uM(1) = 1e-6; % doesnt look like the prediction changes below this value as it reaches zero

[y_ii, y_se] = tauToRate(vertcat(a1.boundmle), vertcat(a1.boundmle_se));

func = @(beta,x) beta(1) + ((beta(2)-beta(1)) ./ (1+beta(3)./x));
minResponse=min(y_ii(:,4));
maxResponse=max(y_ii(:,4));
EC50_guess = 1;
[coeffs,r,J]=nlinfit(x_uM(:),y_ii(:,4),func,[minResponse maxResponse EC50_guess]);
CI = nlparci(coeffs,r,'jacobian',J);
se = (coeffs(end) - CI(3,1)) / 1.96;

disp(['EC50: ', num2str(coeffs(3)), ' ± ', num2str(se)])

xrange = logspace(-6, 2, 500);
xrange = xrange(:);
[Ypred,delta] = nlpredci(func, xrange, coeffs, r, J);
Ypred_ci = delta*[-1 1]+Ypred;
fig3 = figure;
hold on
xrange = log10(xrange);
patch([xrange' fliplr(xrange')], [Ypred_ci(:,2)' fliplr(Ypred_ci(:,1)')], [0 0.4470 0.7410], 'FaceColor', [0 0.4470 0.7410], 'EdgeColor','none', 'FaceAlpha', 0.1);
plot(xrange, Ypred, '-', 'Color', [0 0.4470 0.7410])
scatter(log10(x_uM), y_ii(:,4), 15, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w')
errorbar(log10(x_uM), y_ii(:,4), y_se(:,4), y_se(:,4), 'LineStyle', 'None', 'Color', 'k', 'MarkerFaceColor', [0.3010 0.7450 0.9330], 'CapSize', 0, 'Marker', 'none');

plot(log10(x_uM), y_ii(:,3), '--', 'Color', 'k')
scatter(log10(x_uM), y_ii(:,3), 15, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.3010 0.7450 0.9330])
errorbar(log10(x_uM), y_ii(:,3), y_se(:,3), y_se(:,3), 'k', 'CapSize', 0,...
    'Marker', 'none', 'LineStyle', 'none')
ylabel('kapparent (s-1)')
xlabel('[Branaplam] (µM)')
ylim([1e-4, 1e-1])
set(gca,'yscale','log')

% Break axes
set(gca,'tickdir', 'out', 'fontSize', fontSize, 'box', 'off', 'Layer', 'Top', 'XColor', 'k','YColor', 'k',  'TickLength', [0.03, 0.03])
fig3.Units = 'inches';
fig3.Position = [0, 0, 1.75, 1.5];
xlim([-6.25, 1])
xticks(-6:1:-2)

xticklabels({'10^{-1}', '10^{0}', '10^{1}', '10^{2}'})
set(gcf,'color','w');

legend('off')
h = breakxaxis([-5.75, -1.25]);
h.leftAxes.XTick = -6;
h.leftAxes.XTickLabel = '0';
h.rightAxes.XTick = [-1,0,1,2];
print(gcf, [figurePath, '9bp-1A_Branaplam-EC50-koff2_log_break'],'-dpdf')

%%
%% Bound dwell fit with hill 
x_uM = vertcat(a1.field_value);
x_uM(1) = 1e-6; % doesnt look like the prediction changes below this value as it reaches zero
y_ii = vertcat(a1.boundmle);
y_se = vertcat(a1.boundmle_se);


func = @(beta,x) beta(1) + ((beta(2)-beta(1)) ./ (1+beta(3)./x));
minResponse=min(y_ii(:,4));
maxResponse=max(y_ii(:,4));
EC50_guess = 1;
[coeffs,r,J]=nlinfit(x_uM(:),y_ii(:,4),func,[minResponse maxResponse EC50_guess]);
CI = nlparci(coeffs,r,'jacobian',J);
se = (coeffs(end) - CI(3,1)) / 1.96;

disp(['EC50: ', num2str(coeffs(3)), ' ± ', num2str(se)])

xrange = logspace(-6, 2, 500);
xrange = xrange(:);
[Ypred,delta] = nlpredci(func, xrange, coeffs, r, J);
Ypred_ci = delta*[-1 1]+Ypred;
fig3 = figure;
hold on
xrange = log10(xrange);
patch([xrange' fliplr(xrange')], [Ypred_ci(:,2)' fliplr(Ypred_ci(:,1)')], [0 0.4470 0.7410], 'FaceColor', [0 0.4470 0.7410], 'EdgeColor','none', 'FaceAlpha', 0.1);
plot(xrange, Ypred, '-', 'Color', [0 0.4470 0.7410])
scatter(log10(x_uM), y_ii(:,4), 15, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w')
errorbar(log10(x_uM), y_ii(:,4), y_se(:,4), y_se(:,4), 'LineStyle', 'None', 'Color', 'k', 'MarkerFaceColor', [0.3010 0.7450 0.9330], 'CapSize', 0, 'Marker', 'none');

plot(log10(x_uM), y_ii(:,3), '--', 'Color', 'k')
scatter(log10(x_uM), y_ii(:,3), 15, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.3010 0.7450 0.9330])
errorbar(log10(x_uM), y_ii(:,3), y_se(:,3), y_se(:,3), 'k', 'CapSize', 0,...
    'Marker', 'none', 'LineStyle', 'none')
ylabel('kapparent (s-1)')
xlabel('[Branaplam] (µM)')
%ylim([1e-4, 1e-1])
set(gca,'yscale','log')

% Break axes
set(gca,'tickdir', 'out', 'fontSize', fontSize, 'box', 'off', 'Layer', 'Top', 'XColor', 'k','YColor', 'k',  'TickLength', [0.03, 0.03])
fig3.Units = 'inches';
fig3.Position = [0, 0, 1.75, 1.5];
xlim([-6.25, 1])
xticks(-6:1:-2)

xticklabels({'10^{-1}', '10^{0}', '10^{1}', '10^{2}'})
set(gcf,'color','w');

legend('off')
h = breakxaxis([-5.75, -1.25]);
h.leftAxes.XTick = -6;
h.leftAxes.XTickLabel = '0';
h.rightAxes.XTick = [-1,0,1,2];
print(gcf, [figurePath, '9bp-1A_Branaplam-EC50-tau-bound-2_log_break'],'-dpdf')


%% Event correlation of 1 nM. Bound i to Bound i+1
k = 4; % only using index=4
markerSize = 5;
idx = [2, 4, 6];

plot_data_J = {}; 
for ii = 1:3
    i = idx(ii);
    x = a1(i).field_value;
    events = {};
    for j = 1:numel(a1(i).events)
        events = [events; a1(i).events{j}];
    end
    [dwellPairs, ~, ~, ~, xLabels, yLabels] = dwellTimeCorrelations(events, a1(i).frameRate_s{1}); % only using index=4
    [corr, prob] = corrcoef(log10(dwellPairs{k}));
    plot_data_J{ii} = [log10(dwellPairs{k})];
    nn = numel(dwellPairs{k});
    rr = corr(1,2);
    pp = prob(1,2);
    figName = ['9bp-1A_dwellCorr-BB_', '_', num2str(a1(i).field_value), '-µM'];
    h = plotDwellTimeCorrelations(dwellPairs{k}, 0.3, 0, 4, 1, 0, 1);
    hold on
    if  k == 4
        plot(log10(a1(i).boundmle(3)), log10(a1(i).boundmle(3)), '+k', 'MarkerFaceColor','none', 'MarkerSize', markerSize)
        plot(log10(a1(i).boundmle(4)), log10(a1(i).boundmle(4)), '*k', 'MarkerFaceColor','k', 'MarkerSize', markerSize)
    end
    
    xlabel(xLabels{k})
    ylabel(yLabels{k})
    pbaspect([1,1,1])
    xlim([0,4]);
    ylim([0,4]);
    xticklabels({'10^{0}','10^{1}','10^{2}','10^{3}', '10^{4}'})
    yticklabels({'10^{0}','10^{1}','10^{2}','10^{3}', '10^{4}'})
    % print output
    formatspec  = 'r = %0.2f (N = %i, p = %0.2e)';
    str = sprintf(formatspec, rr, nn, pp);
    disp(str)
    
    % save
    publishFigure(h, ...
        'figureName', [figurePath, figName], ...
        'legendFontSize', 5, ...
        'fontSize', 7, ...
        'tickLength', 0.03,...
        'figureSize', [2, 2]);
end

%%

k =1 ; % only using index=4
markerSize = 5;
idx = [4];
for ii = 1:length(idx)
    i = idx(ii);
    x = a1(i).field_value;
    events = {};
    for j = 1:numel(a1(i).events)
        events = [events; a1(i).events{j}];
    end
    [dwellPairs, ~, ~, ~, xLabels, yLabels] = dwellTimeCorrelations(events, a1(i).frameRate_s{1}); % only using index=4
    [corr, prob] = corrcoef(log10(dwellPairs{k}));
    nn = numel(dwellPairs{k});
    rr = corr(1,2);
    pp = prob(1,2);
    figName = ['9bp-1A_dwellCorr-BB_', '_', num2str(a1(i).field_value), '-µM'];
    % h = plotDwellTimeCorrelations(dwellPairs{k}, 0.3, 0, 4, 1, 0, 1);
    h = plotDwellTimeCorrelations(dwellPairs{k}, 0.3, 0, 4, 1, 0, 1);
    hold on
    if  k == 4
        plot(log10(a1(i).boundmle(3)), log10(a1(i).boundmle(3)), '+k', 'MarkerFaceColor','none', 'MarkerSize', markerSize)
        plot(log10(a1(i).boundmle(4)), log10(a1(i).boundmle(4)), '*k', 'MarkerFaceColor','k', 'MarkerSize', markerSize)
    end
    xlabel(xLabels{k})
    ylabel(yLabels{k})
    pbaspect([1,1,1])
    xlim([0,4]);
    ylim([0,4]);
    xticklabels({'10^{0}','10^{1}','10^{2}','10^{3}', '10^{4}'})
    yticklabels({'10^{0}','10^{1}','10^{2}','10^{3}', '10^{4}'})
    % print output
    formatspec  = 'r = %0.2f (N = %i, p = %0.2e)';
    str = sprintf(formatspec, rr, nn, pp);
    disp(str)
    
    % save
    publishFigure(h, ...
        'figureName', [figurePath, figName], ...
        'legendFontSize', 5, ...
        'fontSize', 7, ...
        'tickLength', 0.03,...
        'figureSize', [2, 2]);
end



%% CDF at 1 uM Branaplam
K = 4;
[h, plot_data_H] = plotDwellsCDF(a1(K).dwells{2}, a1(K).dwellsmle{2}, 'markerSize', 3);
ylim([0,1])
set(gca,'xscale', 'log')
xlim([1 1e4])
xticks([1, 1e1, 1e2, 1e3, 1e4])
xticklabels({'10^0','10^1','10^2','10^3', '10^4'})
ylabel('Cumulative Probability')
xlabel('Bound Time (s)')
legend({['N = ' num2str(numel(a1(K).dwells{2}))], 'monoexp', 'biexp'}, 'location', 'northwest')
publishFigure(h, ...
    'figureName', [figurePath, '1uM_bound_CDF_MLE'], ...
    'legendFontSize', 6, ...
    'fontSize', 7, ...
    'tickLength', 0.03,...
    'figureSize', [1.75, 1.5]);


% 
% %% Load in -U1C 1- uM data
% idx1 = find(vertcat(analysis.cpd_uM) == 10);
% bootstrap = 0;
% tauGuess{1} = [];
% tauGuess{2} = [0.1, 100, 1000];
% useFitLimits = 1;
% a2 = mergeDataByField(analysis(idx1), 'U1C_nM', 'tauGuess', tauGuess, 'fitLimits', useFitLimits, 'bootstrap', bootstrap);
% N = numel(a2);
% 
% %% merge data into 4 categories 
% % its messy but hey it works. 
% A00 = [analysis1; analysis2]; 
% idx1 = find(vertcat(A00.cpd_uM) == 0);
% idx2 = find(vertcat(A00.rna_nM) == 1);
% idx3 = find(vertcat(A00.cpd_uM) == 10);
% 
% B1 = A00(intersect(idx1, idx2));
% b1 = find(vertcat(B1.U1C_nM) == 0);
% b2 = find(vertcat(B1.U1C_nM) == 100);
% 
% a1 = mergeDataByField(B1(b1), 'U1C_nM', 'tauGuess', {[], [0.5, 50, 100]}, 'fitLimits', 1); 
% a2 = mergeDataByField(B1(b2), 'U1C_nM', 'tauGuess', {[], [0.5, 50, 300]}, 'fitLimits', 1); 
% 
% B2 = A00(intersect(idx2, idx3));
% b3 = find(vertcat(B2.U1C_nM) == 0);
% b4 = find(vertcat(B2.U1C_nM) == 100);
% 
% a3 = mergeDataByField(B2(b3), 'U1C_nM', 'tauGuess', {[], [0.5, 50, 50]}, 'fitLimits', 1); 
% a4 = mergeDataByField(B2(b4), 'U1C_nM', 'tauGuess', {[], [0.5, 50, 1000]}, 'fitLimits', 1); 
% 
% % merge 
% A00 = [];
% A00 = [a1; a3; a2; a4];
% 
% %% Plot violin plot of all four options
% str = cell(4,1);
% X = [];
% Y = [];
% for i = 1:4
%     x = A00(i).dwells{2};
%     y = x*0 + i;
%     X = [X;x];
%     Y = [Y;y];
%     %str{i} = [sprintf('%s%i\n%0.1f', 'N = ', numel(x), A00(i).dwellsmle{2}.monoExpTau), ' s'];
%     str{i} = [num2str(numel(x))];
% end
% 
% h=figure;
% vp = violinplot(log10(X), Y, 'ShowData', false, 'EdgeColor', [0,0,0],...
%     'MedianColor', [1,1,1], 'BoxColor',  [0,0,0], 'ViolinAlpha', 0.75);
% hold on
% % yline(log10(mean( A00(i).dwells{2})),'-r')
% ylim([0,4])
% yticks(0:1:5)
% yticklabels({'10^{0}','10^{1}','10^{2}','10^{3}', '10^{4}'})
% 
% ylabel('Bound Time (s)')
% figname = 'bound_violin_';
% 
% for i = 1:4
%     text(i,4, str{i}, 'FontName', 'Arial','HorizontalAlignment', 'center', 'FontSize', 5);
% end
% % xlabel('[U1C] (nM)')
% xlim([0.5, 4.5])
% publishFigure(h, ...
%     'figureName', [figurePath, figname], ...
%     'legendFontSize', 6, ...
%     'fontSize', 7, ...
%     'tickLength', 0.03,...
%     'figureSize', [2, 1.65]);
% 
% 
% %% Plot overlaid CDF of all conditions
% str = cell(4,1);
% X = [];
% Y = [];
% h = figure; 
% names = {'-U1C, +DMSO', '-U1C, +10µM Branaplam', '+U1C, +DMSO', '+U1C, +10µM Branaplam'};
% for i = 1:4
%     x = A00(i).dwells{2};
%     y = x*0 + i;
%     X = [X;x];
%     Y = [Y;y];
%     % str{i} = [sprintf('%s%i\n%0.1f', 'N = ', numel(x), A00(i).dwellsmle{2}.monoExpTau), ' s'];
%     str{i} = [names{i}, ' N=', num2str(numel(x))]
%     
%     % 
%     
%     [time_s, ~, count] = countCDF(x(:));
%     hold on
%     plot(time_s, count)
%     
% end
% ylim([0,1])
% set(gca,'xscale', 'log')
% xlim([1 1e4])
% xticks([1, 1e1, 1e2, 1e3, 1e4])
% xticklabels({'10^0','10^1','10^2','10^3', '10^4'})
% ylabel('Cumulative Probability')
% xlabel('Bound Time (s)')
% %legend(str, 'location', 'northwest')
% publishFigure(h, ...
%     'figureName', [figurePath, 'effect_of_U1C-Bran_CDF'], ...
%     'legendFontSize', 6, ...
%     'fontSize', 7, ...
%     'tickLength', 0.03,...
%     'figureSize', [1.75, 1.25]);
