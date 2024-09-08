%% Make Figure 03. Effect of U1-C
% 2024-01-31
% David S. White 

clear
close all 

%% set file path 
figurePath = '/Users/work/Documents/Research/White_2023/figures-raw/figure04/'; 

%% genreal 
fontSize = 7;

%% Load data
load('analysis_9bp-1A.mat')

nFiles = numel(analysis);

%% Grab the data at 1 nM 9bp-1A with variable U1-C
idx = zeros(nFiles,1);
for i = 1:nFiles
    if analysis(i).rna_nM == 1 && analysis(i).cpd_uM == 0
        idx(i) = 1;
    end
end
% index data
analysis1 = analysis(idx==1);
% analysis1(11) = [];

% sort data by U1-C_nM
[~, sorted_idx] = sort(vertcat(analysis1.U1C_nM));
analysis1 = analysis1(sorted_idx);

% combine seperate experiments under same conditions (ie, same RNA_nM)
a1 = mergeDataByField(analysis1, 'U1C_nM','fitLimits', 1);
N = length(a1);
X = vertcat(a1.field_value);

%% EC50 

N = length(analysis1);
X = vertcat(analysis1.U1C_nM); 
X(X==0) = 1e-6; 
fractionbound = zeros(N,1);
for i = 1:N
    idx = analysis1(i).nevents > 1; 
    fractionbound(i) = mean(analysis1(i).fractionbound(idx==1));
end


plot_A_raw = [X, fractionbound(:)];

func = @(beta,x) beta(1) + ((beta(2)-beta(1)) ./ (1+beta(3)./x));
minResponse=min(fractionbound(:,1));
maxResponse=max(fractionbound(:,1));
EC50_guess = 1;
[coeffs,r,J]=nlinfit(X(:),fractionbound(:,1),func,[minResponse maxResponse EC50_guess]);
CI = nlparci(coeffs,r,'jacobian',J);
se = (coeffs(end) - CI(3,1)) / 1.96;

disp(['EC50: ', num2str(coeffs(3)), ' ± ', num2str(se)])

xrange = logspace(-6,2, 500);
xrange = xrange(:);
[Ypred,delta] = nlpredci(func, xrange, coeffs, r, J);
Ypred_ci = delta*[-1 1]+Ypred;
fig3 = figure;
hold on
xrange = log10(xrange);
patch([xrange' fliplr(xrange')], [Ypred_ci(:,2)' fliplr(Ypred_ci(:,1)')], 'r', 'FaceColor', 'r', 'EdgeColor','none', 'FaceAlpha', 0.1);
plot(xrange, Ypred, '-r')
% mean fraction bound, error = standard error
x = X;
x(1) = 1e-6; % doesnt look like the prediction changes below this value as it reaches zero
scatter(log10(x), fractionbound, 15, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w')
ylabel('Fraction Bound')
xlabel('[U1-C] (nM)')
ylim([0, 0.6])

plot_A_fit = [xrange(:), Ypred(:), Ypred_ci];

% Break axes
set(gca,'tickdir', 'out', 'fontSize', fontSize, 'box', 'off', 'Layer', 'Top', 'XColor', 'k','YColor', 'k',  'TickLength', [0.03, 0.03])
fig3.Units = 'inches';
fig3.Position = [0, 0, 1.75, 1.5];
ylim([0, 0.6])
xlim([-6.25, 2])
xticks(-6:1:-2)

xticklabels({'10^{-1}', '10^{0}', '10^{1}', '10^{2}'})
set(gcf,'color','w');

legend('off')
h = breakxaxis([-5.75, -1.25]);
h.leftAxes.XTick = -6;
h.leftAxes.XTickLabel = '0';
h.rightAxes.XTick = [-1,0,1,2];
print(gcf, [figurePath, '9bp-1A_U1-C-EC50_log_break'],'-dpdf')

%% Dwells- overlay 
U1C_nM = [0, 0.1, 1, 10, 100];
n = length(U1C_nM);
plot_data_B = {}; 
plot_data_D = {}; 
for i=1:2
    if i == 1
        [cs]= colorGradient([0.7,0.7,0.7], [0.8500 0.3250 0.0980], n);
    elseif i == 2
        [cs]= colorGradient([0.7,0.7,0.7], [0 0.4470 0.7410], n);
    end
    disp(i)
    h = figure;
    hold on
    ls = {};
    for j = 1:length(U1C_nM)
        ii = find(vertcat(a1.field_value) == U1C_nM(j));
        %if i == 1
        %    x = a1(ii).timeToFirst;
        %else
        %    x = a1(ii).dwells{i};
        %end
        % disp(num2str(length(x)));
        x = a1(ii).dwells{i};
        disp(num2str(length(x)))
        
        [x_unique,~,y_cdf] = countCDF(x);
        
        ls{j} = [num2str(a1(ii).field_value), ' nM'];
        % scatter(x_unique, y_cdf, 2, 'Marker', 'o', 'MarkerFaceColor', cs(j,:), 'MarkerEdgeColor',  cs(j,:))
        plot(x_unique, y_cdf, 'color', cs(j,:)) 
        
        if i == 1
            plot_data_B{j} = [x_unique(:), y_cdf(:)];
        else
            plot_data_D{j} = [x_unique(:), y_cdf(:)];
        end
    end
    ylabel('Cumulative Probability')
    if i == 1
        xlabel('Unbound Time (s)')
        figureName = [figurePath, 'unboundCDF'];
    else
        xlabel('Bound Time (s)')
        figureName = [figurePath, 'boundCDF'];
    end
    %xlim([0,1800])
    %xticks(0:600:1800)
    set(gca,'xscale','log')
    xlim([1,1e4])
    legend(ls,'Location', 'northwest')
    xticks([1, 10, 100, 1000, 1e4])
    xticklabels({'10^{0}','10^{1}','10^{2}','10^{3}', '10^{4}'})
    publishFigure(h,...
        'figureName', figureName,...
        'fontSize', fontSize,...
        'legendFontSize', fontSize-1,...
        'tickLength', 0.03,...
        'figureSize', [1.75, 1.5]);
end


%% dwells unbound hill TAU
x_nM = vertcat(a1.field_value);
x_nM(1) = 1e-6;
X = x_nM
y_ii = vertcat(a1.unboundmle)
y_se = vertcat(a1.unboundmle_se)
func = @(beta,x) beta(1) + ((beta(2)-beta(1)) ./ (1+beta(3)./x));
minResponse=min(y_ii(:,1));
maxResponse=max(y_ii(:,1));
EC50_guess = 1;
[coeffs,r,J]=nlinfit(X(:),y_ii(:,1),func,[minResponse maxResponse EC50_guess]);
CI = nlparci(coeffs,r,'jacobian',J);
se = (coeffs(end) - CI(3,1)) / 1.96;

disp(['EC50: ', num2str(coeffs(3)), ' ± ', num2str(se)])

xrange = logspace(-6,2, 500);
xrange = xrange(:);
[Ypred,delta] = nlpredci(func, xrange, coeffs, r, J);
Ypred_ci = delta*[-1 1]+Ypred;
fig3 = figure;
hold on
xrange = log10(xrange);
patch([xrange' fliplr(xrange')], [Ypred_ci(:,2)' fliplr(Ypred_ci(:,1)')], [0.8500 0.3250 0.0980], 'FaceColor', [0.8500 0.3250 0.0980], 'EdgeColor','none', 'FaceAlpha', 0.1);
plot(xrange, Ypred, '-', 'Color', [0.8500 0.3250 0.0980])

plot_data_C_fit = [xrange(:), Ypred(:), Ypred_ci];

x = X;
x(1) = 1e-6; % doesnt look like the prediction changes below this value as it reaches zero
scatter(log10(x), y_ii(:,1), 15, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.8500 0.3250 0.0980])
errorbar(log10(x), y_ii(:,1), y_se(:,1), y_se(:,1), 'k', 'MarkerSize', 0.1, 'CapSize', 0,...
    'Marker', 'none', 'LineStyle', 'none')
ylabel('tau unbound (s)')
xlabel('[U1-C] (nM)')
%ylim([1.75e-3, 4.25e-3])

% Break axes
set(gca,'tickdir', 'out', 'fontSize', fontSize, 'box', 'off', 'Layer', 'Top', 'XColor', 'k','YColor', 'k',  'TickLength', [0.03, 0.03])
fig3.Units = 'inches';
fig3.Position = [0, 0, 1.75, 1.5];
xlim([-6.25, 2])
xticks(-6:1:-2)

xticklabels({'10^{-1}', '10^{0}', '10^{1}', '10^{2}'})
set(gcf,'color','w');

legend('off')
h = breakxaxis([-5.75, -1.25]);
h.leftAxes.XTick = -6;
h.leftAxes.XTickLabel = '0';
h.rightAxes.XTick = [-1,0,1,2];
print(gcf, [figurePath, '9bp-1A_U1-C-EC50-tau-unbound_log_break'],'-dpdf')


%% Bound dwells amplitude as hill slope
x_nM = vertcat(a1.field_value);
x_nM(1) = 1e-6;
X = x_nM
% [y_ii, y_se] = tauToRate(vertcat(a1.unboundmle), vertcat(a1.unboundmle_se));
y_ii = vertcat(a1.boundmle);
y_se = vertcat(a1.boundmle_se);

y_ii = 1-y_ii(:,2);
y_se = y_se(:,2);

func = @(beta,x) beta(1) + ((beta(2)-beta(1)) ./ (1+beta(3)./x));
minResponse=min(y_ii(:,1));
maxResponse=max(y_ii(:,1));
EC50_guess = 1;
[coeffs,r,J]=nlinfit(X(:),y_ii(:,1),func,[minResponse maxResponse EC50_guess]);
CI = nlparci(coeffs,r,'jacobian',J);
se = (coeffs(end) - CI(3,1)) / 1.96;

disp(['EC50: ', num2str(coeffs(3)), ' ± ', num2str(se)])

xrange = logspace(-6,2, 500);
xrange = xrange(:);
[Ypred,delta] = nlpredci(func, xrange, coeffs, r, J);
Ypred_ci = delta*[-1 1]+Ypred;
fig3 = figure;
hold on
xrange = log10(xrange);
patch([xrange' fliplr(xrange')], [Ypred_ci(:,2)' fliplr(Ypred_ci(:,1)')], [0 0.4470 0.7410], 'FaceColor', [0 0.4470 0.7410], 'EdgeColor','none', 'FaceAlpha', 0.1);
plot(xrange, Ypred, '-', 'Color', [0 0.4470 0.7410])

x = X;
x(1) = 1e-6; % doesnt look like the prediction changes below this value as it reaches zero
scatter(log10(x), y_ii(:,1), 15, 'o', 'MarkerEdgeColor', 'w', 'MarkerFaceColor', [0 0.4470 0.7410])
errorbar(log10(x), y_ii(:,1), y_se(:,1), y_se(:,1), 'k', 'MarkerSize', 0.1, 'CapSize', 0,...
    'Marker', 'none', 'LineStyle', 'none')
ylabel('Ampltitude of Tau B2')
xlabel('[U1-C] (nM)')

plot_data_F_Amp_raw = [x_nM, y_ii, y_se];
plot_data_F_Amp_fit = [xrange, Ypred, Ypred_ci];

% Break axes
set(gca,'tickdir', 'out', 'fontSize', fontSize, 'box', 'off', 'Layer', 'Top', 'XColor', 'k','YColor', 'k',  'TickLength', [0.03, 0.03])
fig3.Units = 'inches';
fig3.Position = [0, 0, 1.75, 1.5];
xlim([-6.25, 2])
xticks(-6:1:-2)
xticklabels({'10^{-1}', '10^{0}', '10^{1}', '10^{2}'})
set(gcf,'color','w');
ylim([0,1.1])
yticks(0:0.2:1)

legend('off')
h = breakxaxis([-5.75, -1.25]);
h.leftAxes.XTick = -6;
h.leftAxes.XTickLabel = '0';
h.rightAxes.XTick = [-1,0,1,2];
print(gcf, [figurePath, '9bp-1A_U1-C-EC50-Amplitude_log_break'],'-dpdf')


%% Bound dwells, MLE 
x_nM = vertcat(a1.field_value);
x_nM(1)=1e-6;
cs = colorGradient([1,1,1], [0 0.4470 0.7410], 3);
y_ii = vertcat(a1.boundmle);
y_se = vertcat(a1.boundmle_se);

h1 = figure; hold on
plot(log10(x_nM), y_ii(:,4), '-', 'Color', 'k')
scatter(log10(x_nM), y_ii(:,4), 15, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.3010 0.7450 0.9330])
errorbar(log10(x_nM), y_ii(:,4), y_se(:,4), y_se(:,4), 'LineStyle', 'None', 'Color', 'k', 'MarkerFaceColor', [0.3010 0.7450 0.9330], 'MarkerSize', 0.1, 'CapSize', 0, 'Marker', 'o');

plot(log10(x_nM), y_ii(:,3), '-', 'Color', 'k')
scatter(log10(x_nM), y_ii(:,3), 15, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w')
errorbar(log10(x_nM), y_ii(:,3), y_se(:,3), y_se(:,3), 'Color', 'k', 'LineStyle', 'None', 'MarkerFaceColor', 'w', 'MarkerSize', 0.1, 'CapSize', 0, 'Marker', 'o');

legend({'{\it k_{off}^1}','{\it k_{off}^2}'}, 'location','best')


plot_data_F = [x_nM, y_ii(:,3), y_se(:,3),  y_ii(:,4), y_se(:,4)];
    
% previous label but it messes with axes size.. fix in Adobe
% ylabel('{\it k_{off}} (s^{-1})')
ylabel('tau bound (s)')
xlabel('[U1-C] (nM)')
%ylim([1e-4, 1e-1])
set(gca,'yscale','log')
ax = gca;
%ax.YAxis.Exponent = -2;
ylim([1e1, 1e4])
%yticklabels({'0.4', '0.3', '0.2', '0.1'})

set(gca,'tickdir', 'out', 'FontSize', 7, 'FontName', 'Arial',...
    'box', 'off', 'Layer', 'top', 'XColor', 'k', 'YColor', 'k','TickLength', [0.03, 0.03])
h1.Units = 'inches';
h1.Position = [0, 0, 1.75, 1.5];
xlabel('[U1-C] (nM)')
% xlim([-6.25, 2.25])
xlim([-6.25, 2])
xticks(-6:1:-1)
xticklabels({'10^{-1}', '10^{0}', '10^{1}', '10^{2}'})
set(gcf,'color','w');

legend({'{\it k_{off}^1}','{\it k_{off}^2}'}, 'location','best')
print(gcf, [figurePath, '9bp-1A_vary-U1-C_bound-koff_MLE_legend'], '-dpdf')

h = breakxaxis([-5.75, -1.25]);
h.leftAxes.XTick = -6;
h.leftAxes.XTickLabel = '0';
h.rightAxes.XTick = [-1,0,1,2];
legend({'{\it k_{off}^1}','{\it k_{off}^2}'}, 'location','best')
legend('off')
print(gcf, [figurePath, '9bp-1A_vary-U1-C_bound-koff_MLE'], '-dpdf')

%% Bound Tau Amplitue


%% Dwell times at 1 nM U1-C as CDF
K = 4;
[h, plot_data_E] = plotDwellsCDF(a1(K).dwells{2}, a1(K).dwellsmle{2},'markerSize', 5);
ylim([0,1])
%set(gca,'xscale', 'log')
%xlim([1 10000])
ylabel('Cumulative Probability')
xlabel('Bound Time (s)')
set(gca,'xscale','log')
xlim([1,1e4])
legend(ls,'Location', 'northwest')
xticks([1, 10, 100, 1000, 1e4])
xticklabels({'10^{0}','10^{1}','10^{2}','10^{3}', '10^{4}'})
legend({['N = ' num2str(numel(a1(K).dwells{2}))], 'monoexp', 'biexp'}, 'location', 'northwest')
publishFigure(h,...
    'figureName', [figurePath, '9bp-1A_U1-C-1nM_DMSO_unbound-dwells-MLE'],...
    'fontSize', fontSize,...
    'legendFontSize', fontSize-1,...
    'tickLength', 0.03,...
    'figureSize', [1.75, 1.5]);

h.Children(end).XAxis


%% Bound event correlation at 0, 1, 100 nM U1-C
k = 4; % only using index=4
markerSize = 5;
idx = [1, 4, 8];

plot_data_G = {}; 

koff_fits = mean(vertcat(a1.boundmle)); 
for ii = 1:length(idx)
    i = idx(ii);
    x = a1(i).field_value;
    events = {};
    for j = 1:numel(a1(i).events)
        events = [events; a1(i).events{j}];
    end
    [dwellPairs, ~, ~, ~, xLabels, yLabels] = dwellTimeCorrelations(events, a1(i).frameRate_s{1}); % only using index=4
    [corr, prob] = corrcoef(log10(dwellPairs{k}));
    plot_data_G{ii} = log10(dwellPairs{k});
    nn = numel(dwellPairs{k});
    rr = corr(1,2);
    pp = prob(1,2);
    figName = ['9bp-1A_dwellCorr-BB_', '_', num2str(a1(i).field_value), '-nM'];
    
    
    h = plotDwellTimeCorrelations(dwellPairs{k}, 0.25, 0, 3.5, 1, 0, 1);
    hold on
    
    %plot(log10(a1(i).boundmle(3)), log10(a1(i).boundmle(3)), '+k', 'MarkerFaceColor','none', 'MarkerSize', markerSize, 'Linewidth',1)
    %plot(log10(a1(i).boundmle(4)), log10(a1(i).boundmle(4)), '*k', 'MarkerFaceColor','k', 'MarkerSize', markerSize,'Linewidth',1)
    
    %xline(log10(koff_fits(3)), '-', 'color', [0.7,0.7,0.7])
    %xline(log10(koff_fits(4)), '-', 'color', [0.7,0.7,0.7])
    %yline(log10(koff_fits(3)), '-', 'color', [0.7,0.7,0.7])
    %yline(log10(koff_fits(4)), '-', 'color', [0.7,0.7,0.7])
    
    
    %xlabel(xLabels{k})
    %ylabel(yLabels{k})
    xlabel('Bound_i (s)')
    ylabel('Bound_{i+1} (s)')
    pbaspect([1,1,1])
    xlim([0,4]);
    ylim([0,4]);
    xticks(0:4)
    yticks(0:4)
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
        'figureSize', [1.75, 1.75]);
end

