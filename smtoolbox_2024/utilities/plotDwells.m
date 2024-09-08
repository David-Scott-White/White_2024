function [h, plotData] = plotDwells(dwells, dwellFits, varargin)
% David S. White
% dwhite7@wisc.edu

% Plot results of a fitting dwell times. see fitDwells.m

% Updates:
% 2021-02-04 DSW Wrote code from fitExpDist2.m
% 2022-12-06 DSW updated to handle customs bins of uneven bin width

dwells = dwells(:);

plotData = {};
plotData{1} = dwells;

% Set defaults if not provided
bins = 'auto';              % can be number of bins, edges, or string
normalization = 'count';    % count or PDF work currently. see histogram.m
plotType = 'bar';           % bar or scatter (need to add stairs and log hist, & cummulative)
plotStairs = 0; 
stairsColor = 'k';
countError = [];            % boolean. Default 1 for scatter, 0 for bar
fitError = 0;               % boolean. Default 1 for scatter, 0 for bar * not currently functional do not use!
edgeColor = 'k';            % color of edges (and errorbars)
faceColor = [0.9,0.9,0.9];  % face color (and stairs color)
faceAlpha = 1;
showLegend = 1;             % boolean.
dropBinThreshold = 1;       % do not show bins with bin counts 1-threshold
nExp = [1,2];               % control which exponentials to plot
markerSize = 5;
capSize = 0;
plotSize = [];
emptyBins = 1;              % boolean. If true, merges empty bins
legendLocation = 'northeast';        % see legend

h = figure;
% check varargin and set defaults if needed
for i = 1:2:length(varargin)-1
    switch varargin{i}
        case 'bins'
            bins = varargin{i+1};
        case 'normalization'
            normalization = varargin{i+1};
        case 'plotType'
            plotType = varargin{i+1};
        case 'countError'
            countError = varargin{i+1};
        case 'fitError'
            fitError = varargin{i+1};
        case 'edgeColor'
            edgeColor = varargin{i+1};
        case 'faceColor'
            faceColor = varargin{i+1};
        case 'faceAlpha'
            faceAlpha = varargin{i+1};
        case 'showLegend'
            showLegend = varargin{i+1};
        case 'dropBinThreshold'
            dropBinThreshold = varargin{i+1};
        case 'nExp'
            nExp = varargin{i+1};
        case 'markerSize'
            markerSize = varargin{i+1};
        case 'plotSize'
            plotSize = varargin{i+1};
        case 'fontSize'
            plotSize = varargin{i+1};
        case 'emptyBins'
            emptyBins = varargin{i+1};
        case 'legendLocation'
            legendLocation = varargin{i+1};
        case 'capSize'
            capSize = varargin{i+1};
        case 'plotStairs'
            plotStairs = varargin{i+1};
        case 'stairsColor'
            stairsColor = varargin{i+1};
    end
end
if ~isempty(plotSize)
    h.Units = 'inches';
    h.Position = [3,3, plotSize(1), plotSize(2)];
end

% Local functions for exponential distributions.
% should replace writing functions with dwellFits stored equation...
if isempty(dwellFits.limits)
    expPDF1 = @(t,tau1) exppdf(t,tau1);
    expPDF2 = @(t,A1,tau1,tau2) A1.*exppdf(t,tau1) + (1-A1).*exppdf(t,tau2);
else
    tmin = dwellFits.limits(1);
    tmax = dwellFits.limits(2);
    expPDF1 = @(x,tau1) 1/(exp(-tmin/tau1)-exp(-tmax/tau1)) * (exp(-x/tau1)/tau1);
    expPDF2 = @(x, A1, tau1, tau2) (A1/tau1*exp(-x/tau1) + (1-A1)/tau2*exp(-x/tau2))/...
        ((A1*(exp(-tmin/tau1)-exp(-tmax/tau1))) + ((1-A1)*(exp(-tmin/tau2)-exp(-tmax/tau2))));
end

% Determine bins.
if isempty(bins)
    bins = 'auto';
end
% optimal bins, requires knutils-master for now.
if ischar(bins)
    switch bins
        case 'opt'
            bins = calcoptimalhistbins(dwells', 50, []);
            if bins < 3
                bins = 'auto';
            end
    end
end
if ischar(bins)
    [binCounts, binEdges] = histcounts(dwells, 'BinMethod', bins, 'Normalization', normalization);
else
    % assume bins provided are binCenters
    % cannot get code to work for "count" if unequal bin size. force.
    switch normalization
        case 'count'
            if length(unique(diff(bins))) > 1
                disp('Warning in plotDwells: Unequal bin sizes only valid for PDF. Switching plotType')
                plotType = 'scatter';
                normalization = 'PDF';
            end
    end
    [binCounts, binEdges] = histcounts(dwells, bins, 'Normalization', normalization);
end

binWidth = diff(binEdges);
binCenters = binEdges(1:end-1)+binWidth/2;

% remove empty bins 
% need a variable to control this. 
binCentersMax0 = binCenters(end);
while binCounts(end) == 0
    binCounts(end) = []; 
    binCenters(end) = []; 
    binEdges(end) = [];
    binWidth(end) = [];
end

% area adjustment of fits (returned as PDFs). Compute count error
switch normalization
    case 'count'
        binCountSum = sum(binCounts);
        binError = sqrt(binCountSum*binCounts./binCountSum .* (1-binCounts./binCountSum));
        areaAdjust = sum(binWidth.*binCounts);
        if binCountSum ~= length(dwells)
            disp('Warning in plotDwells: not all dwells included in plot due to bins')
        end
        
    case 'PDF'
        areaAdjust = 1;
        tempCounts = histcounts(dwells, binEdges);
        N = sum(tempCounts);
        Phat = tempCounts./N;
        binError = sqrt((N.*Phat).*(1-Phat))./(binWidth*N);
        % Previous, delete in future. Does not work if bin widths are not
        % uniform
        % binError = sqrt(N*binCounts./sum(binCounts) .* (1-binCounts./sum(binCounts)))./(binWidth.*N);
        if N ~= length(dwells)
            disp('Warning in plotDwells: not all dwells included in plot due to bins');
        end
end

% reformat to columns
binCenters(:);
binCounts(:);
binError(:);

% determine fitBins (more bins for more resolution of the line)
% xValues = linspace(binCenters(1)-binWidth(1)/2,binCenters(end)*1.1, 1000);
xValues = linspace(0,binCentersMax0*2, 1000);
xValues = xValues(:);

% Evaluate all fits % need option to only plot biexp
yValues = xValues*0;
if isfield(dwellFits, 'monoExpTau') && sum(nExp == 1)
    yValues = expPDF1(xValues,dwellFits.monoExpTau).*areaAdjust;
    p = 1;
end

if isfield(dwellFits, 'biExpTau') && sum(nExp == 2)
    yValues = [yValues, expPDF2(xValues, dwellFits.biExpAmp(1), dwellFits.biExpTau(1), dwellFits.biExpTau(2)).*areaAdjust];
    p = 2;
end

% Plot bin error? standard deviation of binomial distribution
if isempty(countError)
    switch plotType
        case 'bar'
            countError = 0;
        case 'scatter'
            countError = 1;
    end
end

% Plot fitError? 95% confidence interval of parameters from fitDwells.m
if isempty(fitError)
    switch plotType
        case 'bar'
            fitError = 0;
        case 'scatter'
            fitError = 1;
    end
end

% Plot Data (currently scatter and bar)

% Merge empty bins? Per AGATHA (Kaur et al, Methods 2019) Empty bins are
% merged from left to right (still not great... )
if ~emptyBins
    while sum(binCounts == 0)
        bin1 = find(binCounts==0,1);
        bin2 = bin1+1;
        newCenter = mean([binCenters(bin1), binCenters(bin2)]);
        newEdges = mean([binEdges(bin1), binEdges(bin2)]);
        binCenters(bin2) = newCenter;
        binEdges(bin2) = newEdges;
        binEdges(bin1) = [];
        binCenters(bin1) = [];
        binCounts(bin1) = [];
        binError(bin1) = [];
    end
end


switch plotType
    
    case 'bar'
        hb =  bar(binCenters, binCounts, 1);
        hb.EdgeColor = edgeColor;
        hb.FaceColor = faceColor;
        hb.FaceAlpha = faceAlpha;
        % remove pesky astericks
        ha = findobj(gca,'Type','line');
        set(ha,'Marker','none');
        if countError
            hold on
            errorbar(binCenters, binCounts, binError, edgeColor, 'linestyle','none',  'HandleVisibility','off')
        end
        if plotStairs
            hold on
            stairs(binCenters-binCenters(1), binCounts, 'color', stairsColor, 'HandleVisibility','off')
        end
        
    case 'scatter'
        % computer error
        if countError
            errorbar(binCenters, binCounts, binError, binError, 'Color', edgeColor,'marker', 'o',...
                'MarkerFaceColor', faceColor, 'MarkerEdgeColor', edgeColor, 'linestyle', 'none',...
                'markerSize', markerSize, 'CapSize', capSize);
            set(gca,'yScale','log');
        else
            plot(binCenters, binCounts,  'Color', edgeColor,'marker', 'o',...
                'MarkerFaceColor', faceColor, 'MarkerEdgeColor', edgeColor, 'linestyle', 'none',...
                'markerSize', markerSize);
            set(gca,'yScale','log');
        end
end
if countError
    plotData{2} = [binCenters', binCounts', binError'];
else
    plotData{2} = [binCenters', binCounts'];
end

hold on;
% plot fits (in reverse order)
if p == 1
    if fitError
        yLower1 = expPDF1(xValues,dwellFits.monoExpCI(1)).*areaAdjust;
        yUpper1 = expPDF1(xValues,dwellFits.monoExpCI(2)).*areaAdjust;
        errorPlot(xValues,yValues(:,1), [yLower1, yUpper1], 'plotColor', [0, 0.4470, 0.7410]);
    else
        plot(xValues,yValues(:,1),'-','LineWidth',1, 'color', [0, 0.4470, 0.7410]);
    end
    
elseif p == 2
    if fitError
        yLower1 = expPDF1(xValues,dwellFits.monoExpCI(1)).*areaAdjust;
        yUpper1 = expPDF1(xValues,dwellFits.monoExpCI(2)).*areaAdjust;
        yLower2 = expPDF2(xValues, dwellFits.biExpAmpCI(1,1),dwellFits.biExpTauCI(1,1), dwellFits.biExpTauCI(1,2)).*areaAdjust;
        yUpper2 = expPDF2(xValues, dwellFits.biExpAmpCI(2,1),dwellFits.biExpTauCI(2,1), dwellFits.biExpTauCI(2,2)).*areaAdjust;
        errorPlot(xValues, yValues(:,1), [yLower1, yUpper1], 'plotColor', [0, 0.4470, 0.7410]);
        errorPlot(xValues, yValues(:,2), [yLower2, yUpper2], 'plotColor', 'r');
        set(gca ,'Layer', 'Top');
        
    else
        if sum(yValues(:,1))> 0
            plot(xValues,yValues(:,1),'--','LineWidth', 1, 'color',[0, 0.4470, 0.7410]); hold on;
        end
        if sum(yValues(:,2))> 0
            plot(xValues,yValues(:,2),'-','LineWidth', 1, 'color',[0.8500 0.3250 0.0980]);
        end
        set(gca ,'Layer', 'Top');
    end
end
plotData{3} = [xValues, yValues];

% drop bin threshold? Currently, modify xlim based on threshold
if dropBinThreshold < 1
    switch normalization
        case 'PDF'
            binCountsTemp = histcounts(dwells, binEdges);
        case 'count'
            binCountsTemp = binCounts;
    end
    dropIdx = find(cumsum(binCountsTemp) ./ sum(binCountsTemp) > dropBinThreshold);
    xlim([0, binCenters(dropIdx(1))+binWidth/2]);
end


% add figure legend
% **** Legend is messed up due to errorBarPlot... Need to correct
if showLegend
    if p==1
        legend(['N = ', num2str(length(dwells))],...
            ['\tau = ', num2str(round(dwellFits.monoExpTau, 1)), ' ± ', num2str(round(dwellFits.monoExpSE, 1)), ' s'], 'Location', legendLocation)
    elseif p == 2
        if sum(yValues(:,1))> 0
            nd = 1; 
            leg = legend(['N = ', num2str(length(dwells))],...
                ['\tau_0 = ', num2str(round(dwellFits.monoExpTau,nd)), ' ± ', num2str(round(dwellFits.monoExpSE,nd)), ' s'],...
                ['A_1 = ', num2str(round(dwellFits.biExpAmp(1),2)), ' ± ', num2str(round(dwellFits.biExpAmpSE(1),2)), ...
                newline, '\tau_1 = ', num2str(round(dwellFits.biExpTau(1),nd)), ' ± ', num2str(round(dwellFits.biExpTauSE(1),nd)), ' s', ...
                newline, '\tau_2 = ', num2str(round(dwellFits.biExpTau(2),nd)), ' ± ', num2str(round(dwellFits.biExpTauSE(2),nd)), ' s'],...
                'Location', legendLocation);
            leg.ItemTokenSize = leg.ItemTokenSize/3;
        else
            nd = 1;
            leg = legend(['N = ', num2str(length(dwells))],...
                ['A_1 = ', num2str(round(dwellFits.biExpAmp(1),2)), ' ± ', num2str(round(dwellFits.biExpAmpSE(1),2)), ...
                newline, '\tau_1 = ', num2str(round(dwellFits.biExpTau(1),nd)), ' ± ', num2str(round(dwellFits.biExpTauSE(1),nd)), ' s', ...
                newline, '\tau_2 = ', num2str(round(dwellFits.biExpTau(2),nd)), ' ± ', num2str(round(dwellFits.biExpTauSE(2),nd)), ' s'],...
                'Location', legendLocation);
            leg.ItemTokenSize = leg.ItemTokenSize/3;
        end
    else
        
    end
end

% Label axes
if countError
    yMin = min(binCounts(binCounts>0))/10;
else
    yMin = 0;
end
switch normalization
    case 'count'
        ylabel('Count')
        yMin = 0;
        yMax = ceil(max(binCounts)+max(binCounts)*0.1);
        ylim([yMin, yMax]);
        yt = yticks; 
        yMax = yt(end)+yt(2);
        
    case 'PDF'
        ylabel('PDF')
        nonZeroIdx = find(binCounts>0);
        ymin = abs(min(binCounts(nonZeroIdx) - binError(nonZeroIdx)));
        n = floor(log10(ymin));
        yMin = round(ymin, -1*n)/10;
        ymax = max(binCounts + binError);
        switch plotType
            case 'scatter'
                if ymax > 1e-1
                    yMax = 1;
                else
                    yMax = 0.1;
                end
                yticks([1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1]);
            case 'bar'
                yMax = ymax*1.1;
        end
end
% ylim([yMin, yMax]);
xlabel('Time (s)');
xlim([0, binCenters(end)*1.1])
set(gca, 'fontname', 'Arial')

end
