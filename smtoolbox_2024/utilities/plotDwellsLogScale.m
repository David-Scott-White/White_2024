function plotDwellsLogScale(t, dtf,bw)

x = log(t);
x0 = log(mean(t));
xValues = linspace(min(x), max(x), 1000); 
xValues=xValues(:);

% pdf 
gx = @(x, x0, a) a*(exp(x-x0 - exp(x-x0)));

% Plot best fit example
bestFit = dtf.LLR(3); 
yy = zeros(length(xValues), bestFit);
if bestFit == 1
    amp = 1; 
    x0 = log(mean(t));
    yy(:,1) = gx(xValues, x0, amp);
else
    for i= 1:2
        amp = dtf.biExpAmp(i); 
        tau = dtf.biExpTau(i);
        x0 = log(tau);
        yy(:,i) = gx(xValues, x0, amp);
    end
end

% Plotting 
h = figure; hold on
[counts, edges] = histcounts(x, 'binWidth', bw); 
binwidth = diff(edges);

% if pdf 
% counts = counts./(binwidth*sum(counts));
bincenters = edges(1:end-1)+binwidth/2;
hb = bar((bincenters), counts, 'hist');

% hb = bar(exp(edges(1:end-1)), counts,'hist')
hb.EdgeColor = 'k';
hb.FaceColor = [0.9,0.9,0.9];
ha = findobj(gca,'Type','line');
set(ha,'Marker','none');

% if pdf 
% yyScaled = yy;
yyScaled = yy*(sum(binwidth.*counts));
if size(yyScaled,2) > 1
    for i = 1:2
        plot((xValues), yyScaled(:,i), '--k', 'linewidth',0.5);
    end
    plot((xValues), sum(yyScaled,2), '-k', 'linewidth',1);
else
    plot((xValues), yyScaled, '-k', 'linewidth',1);
end
%set(gca,'xscale','log', 'tickdir','out')
ylabel('Count');
% xlabel('log(Time (s))')
xticks = h.CurrentAxes.XTick;
xticklabels = h.CurrentAxes.XTickLabels;
for i = 1:length(xticklabels)
    h.CurrentAxes.XTickLabels{i} = sprintf('%0.0f', exp(xticks(i)));
end
xlabel('Time (s)')
