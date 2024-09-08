function newSequence = sampleSequence(sequence, frameRate_s, exposure_s, interval_s)
% -------------------------------------------------------------------------
% Resample a given time series (high temporal resolution) at a given
% frameRate and exposure (ms)
% -------------------------------------------------------------------------

newSequence = [];
if exposure_s > interval_s
    disp('Error: Exposure must be less than or equal to interval')
    return
end

if ~exist('interval_s', 'var') || isempty(interval_s)
    interval_s = exposure_s;
end

totalTime_s = length(sequence)*frameRate_s;
time_s = frameRate_s:frameRate_s:totalTime_s;
time_s = double(time_s(:));

% check for perfect division of numbers...


newSequence = zeros(totalTime_s/interval_s, 1);
timePassed_s = 0;
frameStart = 1;
frameStop = 2;
for i = 1:length(newSequence)
    
    elapsedTime_s = 0 ;
    while round(elapsedTime_s,5) < round(exposure_s,5)
        frameStop  = frameStop + 1;
        elapsedTime_s = time_s(frameStop)-timePassed_s;
    end
    newSequence(i) = sum(sequence(frameStart:frameStop));
    
    if interval_s > exposure_s
        while round(elapsedTime_s,5) < round(interval_s,5)
            frameStop  = frameStop + 1;
            elapsedTime_s = time_s(frameStop)-timePassed_s;
        end
    end
    
    timePassed_s = timePassed_s + elapsedTime_s;
    frameStart = frameStop+1;
    frameStop = frameStart+1;
end
end