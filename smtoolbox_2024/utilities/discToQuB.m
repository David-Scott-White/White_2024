function allEvents = discToQuB(data, samplingRate_ms, dropFirstEvent, dropLastEvent)
% -------------------------------------------------------------------------
% Export DISCO.mat structure to QuB format
% -------------------------------------------------------------------------
%
% Input:
%   data = DISC structure from DISCO
%   samplingRate_ms = numeric. frame rate in milliseconds
%   dropFirstEvent = Boolean. Include/discard first events?
%   dropLastEvent = Boolean. Include/discard last events?
%
% Output:
%   Writes two files for QuB format with user specified name
%
% Requires:
%   scaleTimeSeries.m
%   computeCenters.m
%   findEvents.m
%
% Note, currently assumes 2 observable states (e.g., unbound and bound)
%
% Author: David S. White
% Updated: 2021-10-08
% License: MIT
% -------------------------------------------------------------------------

if nargin < 3
    dropFirstEvent = 1;
    dropLastEvent = 1;
end
if isfield(data, 'rois')
    rois = data.rois;
else
    rois = data;
end
numROIs = size(rois, 1);

allEvents = [];

% create files
[file, path] = uiputfile('*.txt');
file = file(1:end-4);
traceFile = fopen([path, file,'_QuB.txt'], 'w');  % Trace File  (.txt)
idealFile = fopen([path, file,'_QuB.dwt'], 'w');  % Dwell Times (.dwt)


wb = waitbar(0, 'Exporting DISC to Qub...');
segStart = 0;
skipTrace = 0;
for i = 1:numROIs
    stateSequence = rois(i).disc_fit.class;
    events = findEvents(stateSequence); 
    numEvents = size(events, 1);
    
    if numEvents < 2
        skipTrace = 1;
    end
    
    % Drop first event?
    if dropFirstEvent
        if numEvents > 1
            events = events(2:end,:);
            numEvents = numEvents-1;
        else
            skipTrace = 1;
        end
    end
    
    % Drop last event?
    if dropLastEvent
        if numEvents > 1
            events = events(1:end-1,:);
            numEvents = numEvents-1;
        else
            skipTrace = 1;
        end
    end
    
    % last check
    if numEvents < 2
        skipTrace = 1;
    end
    
    if ~skipTrace
        trace = rois(i).time_series;
        trace = trace(events(1,1):events(end,2));
        stateSequence = stateSequence(events(1,1):events(end,2)); 
        if min(stateSequence)==1
            stateSequence = stateSequence - 1;
            events(:,4) = events(:,4)-1;
        end
        
        % scale the intensity trace to the ideal trace (0, 1, 2 au)
        [~,ideal] = computeCenters(trace, stateSequence);
        trace = scaleTimeSeries(trace, ideal);
        
        % write the trace file
        for f = 1:length(trace)
            fprintf(traceFile,'%f\r\n',trace(f)); % write scaled, raw data
        end
        fprintf(traceFile,'\r\n');
        
        % Write ideal file (0 and 1 states)
        fprintf(idealFile,...
            'Segment: %d Dwells: %d Sampling(ms): %d Start(ms): %d ClassCount: 2 %f %f %f %f \r\n',...
            i, numEvents, samplingRate_ms, segStart,[0, 0.01, 1, 0.01]); %  mean value + sd
        frameDurs = events(:,3);
        states = events(:,4);
        allEvents = [allEvents; [states, frameDurs, frameDurs*0+i]]; % states, frameDuration, molecule idx
        for e = 1:numEvents
            fprintf(idealFile,'%d\t%d\r\n', states(e),frameDurs(e)*samplingRate_ms);
            
        end
        fprintf(idealFile,'\r\n');
        segStart = segStart+(sum(frameDurs*samplingRate_ms));
        
    end
    waitbar(i/numROIs, wb);
    
end
close(wb)
end



