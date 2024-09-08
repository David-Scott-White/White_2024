function simulations = simulateEquilibration(kon, koff, startingState, ligandRange_M, numSimulations, time_s, frameRate_s)
% David S. White 
% 2022-11-20
% MIT 

% simulate equilibration time for a given kon and koff assuming data start
% in the off state. 
% only binary for now. 

kd = koff/kon; 
if ~exist('ligandRange_M', 'var') || isempty(ligandRange_M)
    % order of magnitude of kd
    n=floor( log10(kd));
    ligandRange_M = logspace(n-2,n+2, 20);
    
end

if ~exist('startingState', 'var') || isempty(startingState)
    startingState = 0;
end

if ~exist('numSimulations', 'var') || isempty(numSimulations)
    numSimulations = 1000;
end

if ~exist('time_s', 'var') || isempty(time_s)
    %time_s = round(3.5 * 1/koff); 
    time_s =  3600;
end

if ~exist('frameRate_s', 'var') || isempty(frameRate_s)
    frameRate_s = 1;
end


nFrames = time_s / frameRate_s;
numLigandConcentrations = length(ligandRange_M);
simulations = cell(numLigandConcentrations,1); 
fractionBound = zeros(numLigandConcentrations, 1); 

frames = 1:nFrames;
timeToFirstEvent = cell(numLigandConcentrations,1);
fractionBoundTime = zeros(numLigandConcentrations,nFrames);

wb = waitbar(0, 'Simulating...');

% Need to adjust to grab event based on time_s for different frame rates
for i = 1:numLigandConcentrations
    tau_on = 1/koff;
    tau_off = 1/(kon*ligandRange_M(i));
    simulations{i} = zeros(numSimulations,nFrames);
    timeToFirstEvent{i} = zeros(numSimulations,1);
    % on_times = roundToSpecificNumber(exprnd(tau_on, 1e6, 1), frameRate_s);
    % off_times = roundToSpecificNumber(exprnd(tau_off, 1e6, 1), frameRate_s);
    for j = 1:numSimulations
        % first event
        if startingState == 0
            eventTime = roundToSpecificNumber(exprnd(tau_off, 1, 1), frameRate_s);
        else
            eventTime = roundToSpecificNumber(exprnd(tau_on, 1, 1), frameRate_s);
        end
        eventFrame = round(eventTime/frameRate_s);
        if ~isempty(eventFrame)
            if eventTime > 0
                events = [1, eventFrame, eventFrame, startingState];
            else
                events = [0, 0, 0, 0];
                eventFrame = 0;
            end
            timeToFirstEvent{i}(j) = eventTime;
            while events(end,2) <= nFrames
                if events(end,4) == 0
                    eventTime = roundToSpecificNumber(exprnd(tau_on, 1, 1), frameRate_s);
                    eventFrame = round(eventTime/frameRate_s);
                    newEvent(1,1) = events(end,2)+1;
                    newEvent(1,2) = eventFrame + newEvent(1,1);
                    newEvent(1,3) = newEvent(1,2)-newEvent(1,1)+1;
                    newEvent(1,4) = 1;
                else
                    eventTime = roundToSpecificNumber(exprnd(tau_off, 1, 1), frameRate_s);
                    eventFrame = round(eventTime/frameRate_s);
                    newEvent(1,1) = events(end,2)+1;
                    newEvent(1,2) = eventFrame + newEvent(1,1);
                    newEvent(1,3) = newEvent(1,2)-newEvent(1,1)+1;
                    newEvent(1,4) = 0;
                end
                events = [events; newEvent];
            end
            if events(end,1) > nFrames
                events(end,:) = [];
            end
            events(end,2) = nFrames;
            events(end,3) = nFrames - events(end,1);
            if events(1,1) == 0
                events(1,:) = [];
            end
            simulations{i}(j,:) = eventsToSequence(events);
        end
    end
    fractionBoundTime(i,:) = sum(simulations{i})./numSimulations;
    waitbar(i/numLigandConcentrations, wb);
end
close(wb);


%% Plot Equilibration Time as function of ligand concentration 
h1 = figure; hold on 
legendString = cell(numLigandConcentrations,1);
for i = 1:numLigandConcentrations
    plot(frames*frameRate_s, fractionBoundTime(i,:))
    legendString{i} = sprintf('%0.2d',ligandRange_M(i));
end
legend(legendString,'Location', 'bestoutside')
xlabel('Time (s)')
ylabel('Fraction Bound')




