function handles = guiManualAdjustStates(hObject, handles, p)
% p can be -1 or 1
if handles.checkboxAdjustEvents.Value == 1
    i = handles.info.roiIndex;
    j = handles.info.channelIndex;
    if isfield(handles.data.rois(i,j), 'events')
        events = handles.data.rois(i,j).events;
        currentState = events(handles.info.currentEvent, 4);
        minStates = min(events(:,4));
        maxStates = max(events(:,4));
        newState = currentState + p;
        if newState <= maxStates && newState >= minStates
            s1 = events(handles.info.currentEvent, 1);
            s2 = events(handles.info.currentEvent, 2);
            ss = handles.data.rois(i,j).fit.class;
            ss(s1:s2) = newState;
            [comps, ideal, class] = computeCenters(handles.data.rois(i,j).timeSeries, ss);
            handles.data.rois(i,j).fit.components = comps;
            handles.data.rois(i,j).fit.ideal = ideal;
            handles.data.rois(i,j).fit.class = class;
            handles.data.rois(i,j).events = findEvents(class);
            
            if handles.info.currentEvent-1 == 0 
                handles.info.currentEvent = 1;
            else
                handles.info.currentEvent =  handles.info.currentEvent-1;
            end
                
            %handles.info.currentEvent = handles.info.currentEvent;
            guidata(hObject, handles);
            guiPlotTraces(hObject, handles, 1);
            guiPlotTraces(hObject, handles, 2);
            guiUpdateTraceInfo(hObject, handles);
            
        end
    end
end
