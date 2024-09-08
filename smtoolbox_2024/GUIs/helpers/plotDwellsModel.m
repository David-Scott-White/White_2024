 function model = plotDwellsModel(dwells, model, overlayType)
% -------------------------------------------------------------------------
% Overlay simulation of model with observed dwell times
% -------------------------------------------------------------------------
% 
% IN DEVELOPMENT **********************************************************
% 
% input:
%   model = struct;
%       model.Q
%       model.states
%       model.numTraces
%       model.frameRate_s
%       model.totalTime_s
%       model.dwells (will simulated if needed)
% 
%       note-max model is biexponential kinetics
% 
%   dwells = cell(k,1); *WILL NEED TO ADD IN PLOT PARMS
% 
% 
% Requires:
%   simulateQ.m
%   plotDwells.m
% 
% Demo code:
% fps = 10;
% q = [0, 0.25, 0; 0.9, 0, 0.5; 0, 0.3, 0];
% ss = simulateQ(q, [0, 1, 1], 1e4, fps);
% events = findEvents(ss);
% dwells = dwellTime(events,fps);
% 
% model = struct;
% model.Q = [0,0.25,0; 0.9, 0, 0.5; 0, 0.3, 0];
% model.states  = [0, 1, 1];
% model.numTraces = 100;
% model.frameRate_s = 0.1;
% model.totalTime_s = 1e3;
% plotDwellsModel(dwells, model, 'fit')
% 
% Author: David S. White
% Updated: 2021-10-28
% License: MIT
% -------------------------------------------------------------------------
%
if ~exist('overlayType', 'var') || isempty(overlayType)
    overlayType = 'fit';
end

Simulate the model if dwells does not exist
states = unique(model.states);
numStates = length(states); % can replace with batch simulate
if ~isfield(model, 'dwells')
    model.dwells = cell(numStates,1);
    wb = waitbar(0, 'Simulating model...');
    for i = 1:model.numTraces
        ss = simulateQ(model.Q, model.states, model.totalTime_s, 1/model.frameRate_s);
        events = findEvents(ss);
        for j = 1:numStates
            idx = events(:,4) == states(j);
            if ~isempty(idx)
                model.dwells{j} = [model.dwells{j}; events(idx,3)*model.frameRate_s];
            end
        end
        waitbar(i/model.numTraces, wb);
    end
    close(wb);
end

model.dwellFit = cell(numStates,1);
for i = 1:2
    if length(model.dwells{i}) < 1000
        bootstrap = 1;
    else
        bootstrap = 0;
    end
    model.dwellFit{i} = fitDwells(model.dwells{i}, 2, bootstrap, []);
end

% plot simulated dwells or simulated dwell fits?
%REALLY NEEDS CLEANING should be able to loop over multiple models to compare

for i = 1:numStates

    % currently ripped from plotDwellsObsSim.m (from HCN paper)
    [counts, binedges] = histcounts(dwells{i}, 'Normalization','PDF');
    binwidth = binedges(2)-binedges(1);
    bincenters = binedges(1:end-1) + binwidth/2;
    [~, xMax1] = min(abs(bincenters - ceil(max(dwells{i}+binwidth))));

    h = figure; hold on
    %bar(bincenters(1:xMax1), counts(1:xMax1), 1,...
    %    'EdgeColor', 'k', 'FaceColor', [0.9,0.9,0.9],'linewidth',0.5);

    scatter(bincenters(1:xMax1), counts(1:xMax1),...
       'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.9,0.9,0.9],'linewidth',0.5);
   set(gca,'yscale','log');

    switch overlayType
        case 'raw'
            [counts_sim, bindedges_sim] = histcounts(model.dwells{i}, binedges, 'Normalization','PDF');
            binwidth_sim = bindedges_sim(2)-bindedges_sim(1);
            bincenters_sim = bindedges_sim(1:end-1) + binwidth_sim/2;
            [~, xMax1_sim] = min(abs(bincenters_sim - ceil(max(model.dwells{i}+binwidth_sim))));
            plot(bincenters_sim(1:xMax1_sim), counts_sim(1:xMax1_sim),...
                'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.9,0.9,0.9],'linewidth',0.5);

        case 'fit'
            xx = bincenters(1)/2:binwidth/10:bincenters(xMax1)*1.1;
            df = model.dwellFit{i};

            if df.LLR(end)==2
                yPDF = df.expPDF2(xx, df.biExpAmp(1), df.biExpTau(1), df.biExpTau(2));
            else
                yPDF = df.expPDF1(xx, df.monoExpTau);
            end

            yhat = yPDF.*binwidth.*sum(counts(1:xMax1));
            plot(xx, yhat, '-', 'color', [0.8500, 0.3250, 0.0980],...
                'linewidth',2); hold on;

            xlabel('Dwell Time (s)');
            ylabel('PDF')
            title(['State: ', num2str(states(i))])

            legendText = {['N = ', num2str(length(dwells{i}))], 'Prediction'};
            legend(legendText);

    end
end

end

