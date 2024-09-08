%% fit SPR curves. association and dissociation 
% 2024-01-10

close all 
clear 

%% Set path
figurePath = '/Users/work/Documents/Research/White_2023/figures-raw/figure02/'; 

%% SPR DATA from REMIX. 
% Load raw data, refit with single exponetial... not a good fit for the
% 11bp... 

% -1A
data_1A = readmatrix('20230922 - SPR Raw Data - 10nM U1snRNP binding to -1A bulge in the presence of Branaplam - ABA Co-inject format.xlsx');

% WT 
data_11 = readmatrix('20230922 - SPR Raw Data - 10nM U1snRNP binding to Match RNA in the presence of Branaplam - ABA Co-inject format.xlsx');

%% Sensorgram example of -1A
bran_uM_plot = [0, 0.15625, 0.3125, 0.625, 1.25, 2.5, 5];
cs = colorGradient([0.7,0.7,0.7], [0.4940, 0.1840, 0.5560], numel(bran_uM_plot));
idx = [8, 18, 14, 10, 34, 38, 44];
h = figure; hold on
ls = cell(numel(bran_uM_plot), 1);
plot_data = []; 
for i = 1:length(idx)
    ii = idx(i);
    plot(data_1A(:,ii-1), data_1A(:,ii),'-', 'color', cs(i,:));
    ls{i} = [num2str(round(bran_uM_plot(i),2)), ' µM'];
    plot_data = [plot_data, data_1A(:,ii-1), data_1A(:,ii)]; 
end
xlabel('Time (s)')
ylabel('Response (RU)')
xlim([-10, 720])
xticks(0:120:720)
ylim([-5, 50])
legend(ls)

publishFigure(h,...
    'figureName', [figurePath, 'branaplam_binding_1A_sensorgram'],...
    'legendFontSize', 5, ...
    'fontSize', 7, ...
    'tickLength', 0.03,...
    'figureSize', [1.75, 1.5]);


%% Fit Koff data
dt1  = 1920;
dt2 = 5120; 
koffWT = []; 
koff1A = [];

% hardcoded from provided excel sheets
plot_data_2 = []; 
bran_uM = [0 0 0 0 0.625 0.625 0.3125 0.3125 0.15625 0.15625 0.078125 0.078125 0.0390625 0.0390625 0.01953125 0.01953125 1.25 1.25 2.5 2.5 5 5];
% WT
for i = 1:2:size(data_11,2)-1
    x = data_11(dt1:dt2,i);
    y = data_11(dt1:dt2,i+1);
    f = fit(x,y, 'exp1');
    koffWT = [koffWT, -f.b];
end

% 1A
for i = 1:2:size(data_1A,2)-1
    x = data_1A(dt1:dt2,i);
    y = data_1A(dt1:dt2,i+1);
    y0 = y./(max(y)-data_1A(102,i+1));
    f = fit(x,y, 'exp1'); 
    koff1A = [koff1A, -f.b];
   
end 

plot_data_2 = [bran_uM(:), koffWT(:), koff1A(:)];

% figure; hold on
% plot(bran_uM, koffWT, '^')
% plot(bran_uM, koff1A, 'o')
% set(gca,'xscale', 'log')
% set(gca,'yscale', 'log')
% ylim([1e-4, 1e-2])
% % ylim([0.1e-3, 2.5e-3] 

%% FIT EC50 of -1A
func = @(beta,x) beta(1) + ((beta(2)-beta(1)) ./ (1+beta(3)./x));
minResponse=min(koff1A(:));
maxResponse=max(koff1A(:));
EC50_guess = 1;
[coeffs1A,r1A,J1A]=nlinfit(bran_uM(:),koff1A(:),func,[minResponse, maxResponse EC50_guess]);
CI1A = nlparci(coeffs1A,r1A,'jacobian',J1A);
se1A = (coeffs1A(end) - CI1A(3,1)) / 1.96;
disp(['EC50: ', num2str(coeffs1A(3)), ' ± ', num2str(se1A)])


%% Plot data
% EC50 Log break
xrange = logspace(-6, 1, 1e3)';
xrange(xrange>5)=[];
xrange(end) = 5;
[Ypred,delta] = nlpredci(func,xrange,coeffs1A,r1A,J1A);
Ypred_ci = delta*[-1 1]+Ypred;
fig3 = figure;
hold on

xrange = log10(xrange);
patch([xrange' fliplr(xrange')], [Ypred_ci(:,2)' fliplr(Ypred_ci(:,1)')], '',...
    'FaceColor','r', 'EdgeColor','none', 'FaceAlpha', 0.25);
plot(xrange, Ypred, '-', 'Color', 'r')
xlabel('[Branaplam] (µM)')
ylabel('koff (s-1)')
bran_uM(bran_uM==0) = 1e-6;
scatter(log10(bran_uM), koff1A, 15, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 1)

fit_data = [xrange(:), Ypred(:), Ypred_ci];

% ploy WT data
scatter(log10(bran_uM), koffWT, 15, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 1)
yticks([0,1,2,3]*1e-3)
ylim([0,3*1e-3])
yticklabels({'0', '1', '2','3'})


set(gca,'tickdir', 'out', 'fontSize', 7,'TickLength', [0.03, 0.03],'FontName', 'Arial')
fig3.Units = 'inches';
fig3.Position = [0, 0, 1.75, 1.5];
xlim([-6.25, 1])
xticks(-6:1:-1)
xticklabels({'0.01', '0.1', '1', '10'})
set(gcf,'color','w');
set(gca, 'Layer', 'top', 'XColor', 'k', 'YColor', 'k', 'TickLength', [0.03, 0.03], 'FontName', 'Arial')
legend({'fit', '-1A', 'WT'}, 'location','best')
print(gcf, [figurePath, 'branaplam_binding_1A_koff_EC50_legend'], '-dpdf')
legend('off')

h = breakxaxis([-5.75, -2]);
h.leftAxes.XTick = -6;
h.leftAxes.XTickLabel = '0';
h.rightAxes.XTick = [-2, -1, 0, 1];

xlabel('[Branaplam] (µM)')
print(gcf, [figurePath, 'branaplam_binding_1A_koff_EC50'],'-dpdf')

