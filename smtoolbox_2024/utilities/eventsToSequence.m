function sequence = eventsToSequence(events)
%% Events to State Sequence 
% David S. White 
% 2020-12-28 

% event values from findEvents.mat returned as a state-sequence
% note, affected by if last events are dropped in findEvents

sequence = zeros(events(end,2),1);
for i = 1:size(events,1)
    sequence(events(i,1):events(i,2),1) = events(i,4);
end