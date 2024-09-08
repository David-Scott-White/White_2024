function time_to_first = timeToFirstEvent(events)
% David S. White
% 2023-08-01

% take cell of events {N,1} and return only the duration (frames) of time
% to first binding. 
% The state of time to first binding is determined by minimum value 

events = cell2mat(events); 
first_frame_idx = find(events(:,1) == 1); 
unbound_state_idx = find(events(:,4) == min(events(:,4))); 
time_to_first_idx = intersect(first_frame_idx, unbound_state_idx); 
time_to_first = events(time_to_first_idx, 3); 

end