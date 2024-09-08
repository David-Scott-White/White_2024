function handles = guiSetInitialParamtersTraceViewer(hObject, handles) 
% Set default parameters for trace view. Change as needed
% (option to load a file?) 
handles.info = struct; 

% general 
handles.info.channelIndex = 1;
handles.info.roiIndex = 1;
handles.info.currentEvent = 1;

% Auto trace selection 
handles.info.filter = struct;
handles.info.filter(1).nStates = [1,5];
handles.info.filter(1).SNRsig = [100, inf];
handles.info.filter(1).SNRbkg = [3,inf];
handles.info.filter(1).nStatesOn = 0;
handles.info.filter(1).SNRsigOn = 0;
handles.info.filter(1).SNRbkgOn = 0;

handles.info.filter(2).nStates = [1,5];
handles.info.filter(2).SNRsig = [100, inf];
handles.info.filter(2).SNRbkg = [3, inf];
handles.info.filter(2).nStatesOn = 0;
handles.info.filter(2).SNRsigOn = 0;
handles.info.filter(2).SNRbkgOn = 0;

% apperance (add in options for plot types in future...)
handles.info.channelColor = zeros(2,3); 
handles.info.channelColor(1,:) = [0 0.4470 0.7410]; % Blue
handles.info.channelColor(2,:) = [0.4660 0.6740 0.1880]; % Green
handles.info.histogramBins = 100; 

handles.info.idealColor = [0,0,0]; 
handles.info.idealWidth = 1.5; 

handles.info.highlightColor = [0.8500 0.3250 0.0980]; % Orange
handles.info.highhightWidth = 3; 
handles.info.backgroundColor = [0.8, 0.8, 0.8]; % light grey
handles.info.backgroundIdealColor = [0.5,0.5,0.5]; % dark grey

handles.info.xLabel  = 'Time (s)';
handles.info.ylabel = 'Fluorescence (au)';
handles.info.grid = 'off'; % on or off
handles.info.box = 'off'; % or on 
handles.info.TickDir = 'out'; % or in

handles.info.idealPlotSpotColor = [1 0 0]; 
handles.info.idealPlotSpotWidth = 1; 
handles.info.plotSpotWidth = 80; % NEED AUTOMATED OPTION

guidata(hObject, handles); 

