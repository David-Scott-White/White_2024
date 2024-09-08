function simulation = batchSimulate(Q, states, observationTime_s,...
    frameRate_s, numSimulations)
% -------------------------------------------------------------------------
% Batch simulation
% -------------------------------------------------------------------------
%
% Wrapper around simulateQ to make my life easier.
%
% inputs.
%   see simulateQ
%
% output
%   traces
%   events
%   dwells
%
% David S. White
% 2021-11-03
% MIT
% -------------------------------------------------------------------------
%
%
uniqueStates = unique(states);
k = length(uniqueStates);
traces = cell(numSimulations,1);
events = cell(numSimulations,1);
dwells = cell(k,1);
frameRate_hz = 1/frameRate_s;
if numSimulations > 10
    wb = waitbar(0, 'Running Simulations...');
end
for i = 1:numSimulations
    [traces{i}, p0] = simulateQ(Q, states, observationTime_s, frameRate_hz);
    events{i} = findEvents(traces{i});
    for j = 1:k
        idx = find(events{i}(:,4)==uniqueStates(j));
        if ~isempty(idx)
            dwells{j} = [dwells{j};  events{i}(idx,3)*frameRate_s];
        end
    end
    if numSimulations > 10
        waitbar(i/numSimulations, wb);
    end
end
if numSimulations > 10
    close(wb);
end

simulation = struct;
model.frame_rate_s  = frameRate_s;
simulation.traces = traces;
simulation.events = events;
simulation.dwells = dwells;
simulation.stateOcc = p0;
