figurePath = '/Users/work/Documents/Research/White_2023/figures-raw/figure05/';


%% Average unbound time (data collected at one concentration
names = {'-1C', '+1C', '+2C', '6p-Exon', '6bp-Intron', 'SMN2', 'FOXM1', 'HTT*', '-1A (1 nM)'}


unbound_u1c  =[677.4	226
170.6	40.9
379.6	145.9
144.9	38.8
205.1	51.9
256.4576	38.1368
167	32.7621
103.58	24.78
513.5205 266.4315]

h = figure; hold on 

markers = {'o', 's', 'd', 'o', 's', 'd', '^', 'v', 'o'};
colors = [
    [0 0.4470 0.7410];
    [0.8500 0.3250 0.0980];
    [0.4660 0.6740 0.1880];
    [0.6350 0.0780 0.1840];
    [0.4940 0.1840 0.5560];
    [0 0.4470 0.7410];
    [0.4940 0.1840 0.5560];
    [0.4660 0.6740 0.1880];
    [0,0,0];
    
    ];
min_c = min(unbound_u1c(:));
max_c = max(unbound_u1c(:));
for i = 1:length(unbound_u1c)
    if i > 5
        scatter(unbound_u1c(i,1), unbound_u1c(i,2), 30, markers{i}, 'MarkerEdgeColor', colors(i,:))
    else
         scatter(unbound_u1c(i,1), unbound_u1c(i,2), 30, markers{i}, 'filled', 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', colors(i,:))
    end
end
plot([1e1, 1e3], [1e1, 1e3], '--k')
xlim([1e1, 1e3]);
ylim([1e1, 1e3]);
xticks([1e1, 1e2, 1e3])
yticks([1e1, 1e2, 1e3])
set(gca,'xscale', 'log')
set(gca,'yscale', 'log')
xlabel('Unbound Time (s)')
ylabel('Unbound Time (s)')
pbaspect([1,1,1])
publishFigure(h,...
    'figureName', [figurePath, 'unbound_U1C_effect'],...
    'fontSize', 7,...
    'legendFontSize', 6,...
    'tickLength', 0.03,...
    'legendTokenSize', [10,6],...
    'figureSize', [2, 1.75]);


%%
bound_u1c  = [50.1	325.8
3	5.8
86	171
6.1	8.8
3.6	105.8
5.4925	23.1289
2.0697	16.0338
3.3511	14.4621
65.2881 309.7892];

h = figure; hold on
min_c = min(bound_u1c(:));
max_c = max(bound_u1c(:));
for i = 1:length(bound_u1c)
     if i > 5
        scatter(bound_u1c(i,1), bound_u1c(i,2), 30, markers{i}, 'MarkerEdgeColor', colors(i,:))
    else
         scatter(bound_u1c(i,1), bound_u1c(i,2), 30, markers{i}, 'filled', 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', colors(i,:))
    end
end
plot([1e0, 1e3], [1e0, 1e3], '--k')
%xlim([min_c, max_c]);
%ylim([min_c, max_c]);
%xlim([1e-4, 1e-1])
%ylim([1e-4, 1e-1])
set(gca,'xscale', 'log')
set(gca,'yscale', 'log')
xlim([1e0, 1e3]);
ylim([1e0, 1e3]);
xticks([1e0, 1e1, 1e2, 1e3])
yticks([1e0, 1e1, 1e2, 1e3])
xlabel('Bound Time (s)')
ylabel('Bound Time (s)')
pbaspect([1,1,1])
publishFigure(h,...
    'figureName', [figurePath, 'bound_U1C_effect'],...
    'fontSize', 7,...
    'legendFontSize', 6,...
    'tickLength', 0.03,...
    'legendTokenSize', [10,6],...
    'figureSize', [2, 1.75]);

legend(names)

%xlabel('{\rangle \langle}')
publishFigure(h,...
    'figureName', [figurePath, 'koff_U1C_effect_legend'],...
     'fontName', 'Helvetica',...
    'fontSize', 7,...
    'legendFontSize', 6,...
    'tickLength', 0.03,...
    'legendTokenSize', [10,6],...
    'figureSize', [2, 1.75]);


%%
figure;
scatter(unbound_u1c(:,2)./unbound_u1c(:,1), bound_u1c(:,2)./bound_u1c(:,1))

%% 
[mdl, RMSE, R2] = linearRegression(unbound_u1c(:,2), unbound_u1c(:,1))
plotLinearRegression(mdl, unbound_u1c(:,2), unbound_u1c(:,1))

%%
