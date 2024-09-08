function allEvents = intervalsToQuB(Intervals, samplingRate_ms, dropFirstEvent, dropLastEvent, file, path)
% -------------------------------------------------------------------------
% Export Intervals.mat structure to QuB format
% -------------------------------------------------------------------------
%
% Input:
%   Intervals = Intervals structure from imscroll
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
%
% Make pop up window??

if nargin < 3
    dropFirstEvent = 1;
    dropLastEvent = 1;
end
if isfield(Intervals, 'AllTracesCellArray')
    traceArray = Intervals.AllTracesCellArray;
else
    traceArray = Intervals;
end
numTraces = size(traceArray, 1);

allEvents = [];

% create files
if ~exist('file', 'var') || isempty(file)
    [file, path] = uiputfile('*.txt');
end
if file
    fileTemp = file(1:end-4);
    if dropFirstEvent
        fileTemp = [fileTemp, '_f0'];
    else
        fileTemp = [fileTemp, '_f1'];
    end
    
    if dropLastEvent
        fileTemp = [fileTemp, '_l0'];
    else
        fileTemp = [fileTemp, '_l1'];
    end
    file = [fileTemp, '.txt'];
    
    if ~file
        disp('intervalsToQuB Terminated: Need Export File');
    end
    
    file = file(1:end-4);
    traceFile = fopen([path, file,'_QuB.txt'], 'w');  % Trace File  (.txt)
    idealFile = fopen([path, file,'_QuB.dwt'], 'w');  % Dwell Times (.dwt)
    
    wb = waitbar(0, 'Exporting Intervals to Qub...');
    segStart = 0;
    
    for i = 1:numTraces
        skipTrace = 0;
        trace = traceArray{i,13}(:,3);
        events = traceArray{i,10};
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
        if numEvents < 1 % changed from 2
            skipTrace = 1;
        else
            skipTrace = 0;
        end
        
        if ~skipTrace
            
            % build idealized trace from events (Intervals format)
            stateSequence = trace*0;
            for j = 1:numEvents
                stateSequence(events(j,2):events(j,3)) = events(j,1);
            end
            frameStart = events(1,2);
            frameStop = events(end,3);
            trace = trace(frameStart:frameStop);
            stateSequence = stateSequence(frameStart:frameStop);
            
            % check for extra states (-2, -3, 2, 3)
            stateSequence(stateSequence ==-2) = 0;
            stateSequence(stateSequence == 2) = 0;
            stateSequence(stateSequence ==-3) = 1;
            stateSequence(stateSequence == 3) = 1;
            
            % scale the intensity trace to the ideal trace (0, 1, 2 au)
            [~,ideal] = computeCenters(trace, stateSequence);
            if numEvents > 1
                trace = scaleTimeSeries(trace, ideal);
            else
                trace = trace./ideal;
                if stateSequence(1) == 0
                    trace = trace - 1;
                end
            end
            events = findEvents(stateSequence); % Bit ugly here for now..
            
            
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
        waitbar(i/numTraces, wb);
        
    end
    close(wb)
else
    disp('Error: Need a file to write')
end



