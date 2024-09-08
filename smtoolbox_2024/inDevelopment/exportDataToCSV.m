%% export analysis to csv
% David S. White
% 2023-08-31

%% intended to export the results of smVideoProcessing &/or smTraceViewer to csv
% creates two seperate dirs, one for each channel (ie RR, GG)
% format 
% file 1 [molecule frame, time_s, fluorescence, class, ideal]; 
% file 2 selected vs non selected [N,1]
% file 3 events 
% [molecule number, start frame, stop frame, frame duration, time duration,
% state class, state mu]; 
% file 4, info? 
% only exports selected data!

%% Load data structure
[file,path] = uigetfile('.mat');
load([path,file]); 

%% Grab selected data
selection_idx  = vertcat(data.rois(:,end).status);
rois = data.rois(find(selection_idx==1), :);
[nrois, nchannels] = size(rois); 

%% Create directory
filepath = [path,file(1:end-4),'_export/'];
if ~exist(filepath, 'dir')
    mkdir(filepath)
end
names = {'channel_1', 'channel_2'}; 

%% Write events (could be cleaner but this is still fast so who cares)
event_headers = {'molecule', 'start_frame', 'stop_frame', 'duration_frames', 'state'};
for i = 1:nchannels
    allEvents = []; 
    frameRate_s = round(mean(diff(data.time{i}(:,end))),3);
    for j = 1:nrois
        events = rois(j,i).events; 
        if ~isempty(events)
            nevents = size(events,1); 
            allEvents = [allEvents; zeros(nevents,1)+j,events];
        else
            allEvents = [allEvents; j, NaN, NaN, NaN, NaN];
        end
    end
    T = array2table(allEvents,'VariableNames',event_headers);
    writetable(T, [filepath, names{i}, '_events.csv']);
    % writematrix(allEvents, [filepath, names{i}, '_events.csv']);
    numStates = unique(allEvents(:,end));
    numStates(isnan(numStates)) = []; 
    % write dwells
    for j = 1:numel(numStates)
        dwells = allEvents(allEvents(:,end)==numStates(j), 4) * frameRate_s; 
        writematrix(dwells, [filepath, names{i}, '_state_',num2str(numStates(j)) '_dwells.csv']);
    end
end

%% write traces (with idealization if exists)
 headers = {'molecule', 'time_s', 'time_series', 'fit_class', 'fit_ideal'};
for i = 1:nchannels
    X = []; 
    time_s = round(data.time{i}(:,end)); % assume time is always at least >= 1 second exposure
    for j = 1:nrois
        timeSeries = rois(j,i).timeSeries;
        nFrames = numel(timeSeries);
        if ~isempty(rois(j,i).fit)
            fit_class =  rois(j,i).fit.class;
            fit_ideal = rois(j,i).fit.ideal;
        else
            fit_class = NaN(nFrames,1);
            fit_ideal = NaN(nFrames,1);
        end
        X = [X; zeros(nFrames,1)+j, time_s(1:nFrames), timeSeries, fit_class, fit_ideal];
    end
    T = array2table(X,'VariableNames',headers);
    writetable(T, [filepath, names{i}, '_traces.csv']);
end

%% 
