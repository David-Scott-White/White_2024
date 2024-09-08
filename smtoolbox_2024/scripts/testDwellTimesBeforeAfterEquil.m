kon = 1.1e6; 
koff = 2.1e-3;
startingState = 0; 
ligandRange_M = 1e-9;
numSimulations = 10000; 
time_s = 3600;
frameRate_s = 5; 
simulations = simulateEquilibration(kon, koff, startingState, ligandRange_M, numSimulations, time_s, frameRate_s);

hold on
equilTime = 3.5*1/(kon*ligandRange_M + koff)
xline(equilTime)

disp(['Expected tauOff: ', num2str(1/(kon*ligandRange_M))])
disp(['Expected tauOn: ', num2str(1/(koff))])

%% Parse each molecule to get events. dwell times of bound, time to first, middle
equilFrame = round(equilTime/frameRate_s)
dwells = localGetDwells(simulations{1}(:,equilFrame:end), frameRate_s);

%% Fit all dwells
dwells_off0_fit = fitDwells(dwells{1}, 1, 0, []);
dwells_off1_fit = fitDwells(dwells{2}, 1, 0, []);
dwells_on_fit = fitDwells(dwells{3}, 1, 0, [frameRate_s, time_s]);

%% Plot fits 
plotDwells(dwells{1},dwells_off0_fit); title('dwells off First');
plotDwells(dwells{2},dwells_off1_fit); title('dwells off Middle');
plotDwells(dwells{3},dwells_on_fit); title('dwells on');

%% local function
function dwells = localGetDwells(simulations, frameRate_s)
% first off; all off after first, bound times
dwells = cell(3,1);
numSimulations = size(simulations,1);
for i = 1:numSimulations
    x = simulations(i,:);
    if sum(x) > 0
        x = x(:);
        events = findEvents(x);
        nevents = size(events,1);
        dwells{1} = [dwells{1}; events(1,3)*frameRate_s];
        for j = 2:nevents
            if events(j,4) == 0
                dwells{2} = [dwells{2}; events(j,3)*frameRate_s];
            else
                dwells{3} = [dwells{3}; events(j,3)*frameRate_s];
            end
        end
    end
end
end
