function [allEvents, out] = exportToQuB(allTraces, allIdeal, samplingRate_ms, file, path)
% -------------------------------------------------------------------------
% Export Data to QuB
% -------------------------------------------------------------------------
% A general function for writing data to QuB
%
% input:
%   allTraces = cell(N,1)
%   allIdeal = cell (N,1) (disc_fit.class)
%   samplingRate_ms = numeric.
%   dropFirstEvent = bool
%   dropLastEvent = bool
%
% David S. White
% 2021-10-21
% MIT
% -------------------------------------------------------------------------

numTraces = size(allTraces, 1);
allEvents = [];

out = [];
if ~exist('file', 'var') || isempty(file)
    [file, path] = uiputfile('*.txt');
end

if ~file
    disp('intervalsToQuB Terminated: Need Export File');
else
    
    if strcmp(file, '.txt')
        file = file(1:end-4);
    end
    traceFile = fopen([path, file,'_QuB.txt'], 'w');  % Trace File  (.txt)
    idealFile = fopen([path, file,'_QuB.dwt'], 'w');  % Dwell Times (.dwt)
    
    wb = waitbar(0, 'Exporting Intervals to Qub...');
    segStart = 0;
    
    for i = 1:numTraces
        trace = allTraces{i};
        ideal = allIdeal{i};
        events = findEvents(ideal); % includes first and last!
        numEvents = size(events, 1);
        traceScaled = scaleTimeSeries(trace, ideal);
        
        temp = struct;
        temp.trace = traceScaled;
        temp.ideal = ideal;
        temp.events = events;
        
        % store output for debug
        if isempty(out)
            out = temp;
        else
            out = [out; temp];
        end
        
        % write the trace file
        for f = 1:length(traceScaled)
            fprintf(traceFile,'%f\r\n',traceScaled(f)); % write scaled, raw data
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
        
        waitbar(i/numTraces, wb);
        
    end
    close(wb)
    fclose(traceFile); 
    fclose(idealFile);
end
