function dwells = dwellTime(events,frameRateHz);
% David S. White
% dwhite7@wisc.edu
% updates: 
% 20-02-06 DSW Wrote the function 

% Summary: 
% dwells = cell of [numStates, 1] of event durations adjusted by
% frameRateHz

% input: 

% output: 

% 
if isempty(events)
    dwells = [];
    return
end
states = unique(events(:,4));
nStates = length(states); 
dwells = cell(nStates,1); 
for n = 1:nStates
    idx = find(events(:,4) == states(n));
    dwells{n} = events(idx,3) / frameRateHz;
end