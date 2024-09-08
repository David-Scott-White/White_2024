%% MST data from REMIX
% file 1. MST of branaplam binding

%% Figure Path

figurePath = '/Users/work/Documents/Research/White_2023/figures-raw/figure02/'; 

X = readmatrix('20231107 - MST Raw Data - Branaplam binding experiments - updated with U1C effects.xlsx');

bran_uM = [X(1:16,9); X(1:16,9)]/1000; 
zeroU1C = [X(1:16,10); X(1:16,11)]; 
plusU1C = [X(1:16,12); X(1:16,13)]; 

% normalize to realtive changes in fluoresence
plusU1C = plusU1C-min(plusU1C);
zeroU1C = zeroU1C-min(zeroU1C);

func = @(beta,x) beta(1) + ((beta(2)-beta(1)) ./ (1+beta(3)./x));
minResponse=min(plusU1C);
maxResponse=max(plusU1C);
EC50_guess = 1;
[coeffs,r,J]=nlinfit(bran_uM,plusU1C ,func,[minResponse, maxResponse EC50_guess]);
CI = nlparci(coeffs,r,'jacobian',J);
se = (coeffs(end) - CI(3,1)) / 1.96;
disp(['KD: ', num2str(coeffs(3)), ' ± ', num2str(se)])

% EC50 Log break
xrange = logspace(-4, 2, 500)';
[Ypred,delta] = nlpredci(func,xrange,coeffs,r,J);
Ypred_ci = delta*[-1 1]+Ypred;
h = figure;
hold on
scatter(bran_uM, zeroU1C, 15, '^', 'MarkerFaceColor', [0.5,0.5,0.5], 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha',1)

plot_data_noU1C = [bran_uM(:), zeroU1C(:)];

% xrange = log10(xrange);
patch([xrange' fliplr(xrange')], [Ypred_ci(:,2)' fliplr(Ypred_ci(:,1)')], '',...
    'FaceColor',[0.4940, 0.1840, 0.5560], 'EdgeColor','none', 'FaceAlpha', 0.25);
plot(xrange, Ypred, '-', 'Color', [0.4940, 0.1840, 0.5560])

plot_data_fit = [xrange(:), Ypred(:), Ypred_ci];

%y1 = mean(U1_1A,2); 
%ysd = std(U1_1A')/sqrt(2);
% plot(x,y,'ok', 'MarkerFaceColor', 'w', 'MarkerSize', 4)
% errorbar(bran_uM, y1, ysd, ysd, 'ok', 'MarkerFaceColor', 'w', 'MarkerSize', 4)
scatter(bran_uM, plusU1C, 15, 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 1)
legend({}, 'Location', 'northwest')
set(gca,'xscale','log')
xlim([min(bran_uM), max(bran_uM)])
ylim([-1, 15]); 

plot_data_plusU1C = [bran_uM(:), plusU1C(:)];

ylabel('\Delta Fluorescence (AU)') 
% need to double check ylabel
xlabel('[Branaplam] (µM)')

xticks([1e-4, 1e-3 ,1e-2, 1e-1, 1e0, 1e1,1e2])

publishFigure(h,...
    'figureName', [figurePath, 'MST_branaplam-binding'],...
    'legendFontSize', 6, ...
    'fontSize', 7, ...
    'tickLength', 0.03,...
    'figureSize', [2, 1.85]);
