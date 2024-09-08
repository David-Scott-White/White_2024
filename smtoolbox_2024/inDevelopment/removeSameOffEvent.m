function remove_idx = removeSameOffEvent(rois)
% David S. White 
% 2023-07-21

% takes rois by smTraceViewer, Assumes 2 channels, where channel 1 is the
% surface AOI; 
% Goal is to find "off" events that are identical in each cahnnel to remove
% (ie molecule fell off surface, so signal is lost in both channels). 


N = length(rois); 
toRemove = zeros(N,1);
for i = 1:N
    if ~isempty(rois(i,1).fit) && ~isempty(rois(i,2).fit) && rois(i,2).status == 1
        events1 = rois(i,1).events;
        events2 = rois(i,2).events;
        if size(events1,1) > 1 && size(events2,1) > 1
            % assume channel 1 is surface molecule
            off_idx = events1(:,4) == min(events1(:,4));
            off_frames = events1(off_idx,1);
            % check if a transition happens in the other channel 
            if ~isempty(intersect(events2(:,1), off_frames))
                toRemove(i) = 1;
            end
        end
    end
end
remove_idx = find(toRemove==1); 
end