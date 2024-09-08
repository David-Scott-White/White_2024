function handles = guiInitSMTraceViewer(hObject, handles)
% call function when data is provided as handles.data


cla(handles.axes1) ;
cla(handles.axes2) ;
cla(handles.axes3) ;
cla(handles.axes4) ;

% set some defaults
handles.info.currentEvent = 1;
handles.info.channelIndex = 1;
handles.info.roiIndex = 1;

% init AOI slider (need keyboard key as well)
[nrois, nchannels] = size(handles.data.rois);
handles.info.nchannels = nchannels;
handles.info.nrois = nrois;
handles.data.idealParm = cell(2,1);

if ~isfield(handles.data.rois, 'status')
    [handles.data.rois.status] = deal(0);
end
handles.info.nselected = [0,0]; 
for i = 1:handles.info.nchannels
    handles.info.nselected(i) = sum(vertcat(handles.data.rois(:,i).status)>0);
end

% Get general info
if ~isfield(handles.data.rois(1,1), 'snr') || ~isfield(handles.data.rois(1,1), 'snb')
    handles = guiEstimateAllSNR(hObject, handles);
end

handles.sliderAOI.Min = 1;
handles.sliderAOI.Max = nrois;
if handles.sliderAOI.Max > 1
    handles.sliderAOI.SliderStep = [1/(handles.sliderAOI.Max-1) , 1/(handles.sliderAOI.Max-1)];
else
    handles.sliderAOI.SliderStep = [0, 0];
end
handles.sliderAOI.Value = 1;
handles.textCurrentAOI.String = num2str(handles.sliderAOI.Value);
handles.textTotalAOI.String = num2str(handles.sliderAOI.Max);

% note, might be in old DISCO format, need to update over time

% update handles and plot
guidata(hObject, handles);
if nchannels == 1
    handles.popupmenuView1.Value = 1;
    handles.popupmenuView2.Value = 2;
    
else
    handles.popupmenuView1.Value = 1;
    handles.popupmenuView2.Value = 3;
end

guiPlotTraces(hObject, handles, 1);
guiPlotTraces(hObject, handles, 2);
guiUpdateTraceInfo(hObject, handles);
