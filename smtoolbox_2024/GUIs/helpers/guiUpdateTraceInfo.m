function handles = guiUpdateTraceInfo(hObject, handles)

% update the trace info on display when trace is idealized or swithched
%
%
% David S. White
% 2022-02-08
% MIT
% -------------------------------------------------------------------------
% get index
% get channel (from top of GUI)
% update

chaIndex = handles.info.channelIndex;
roiIdx = handles.info.roiIndex;

handles.textCurrentChannel.String = num2str(chaIndex);

if isfield(handles, 'data')
    if isfield(handles.data, 'rois')
        roi = handles.data.rois(roiIdx, chaIndex);
        if isfield(roi, 'status')
            if isempty(roi.status)
                roi.status = 0;
            end
            handles.textStatus.String = num2str(roi.status);
        else
            handles.textStatus.String = 'n/a';
        end
        
        if isfield(roi, 'fit') && ~isempty(roi.fit)
            handles.textNumStates.String = num2str(size(roi.fit.components,1));
            handles.checkboxAdjustEvents.Enable = 'on';
        else
            handles.textNumStates.String = 'n/a';
            handles.checkboxAdjustEvents.Value = 0;
             handles.checkboxAdjustEvents.Enable = 'off';
        end
        
        if isfield(roi, 'events')
            handles.textNumEvents.String = num2str(size(roi.events,1));
        else
            handles.textNumEvents.String = 'n/a';
        end
        
        if isfield(roi, 'snr')
            handles.textSignalToNoise.String = num2str(round(roi.snr, 1));
        else
            handles.textSignalToNoise.String = 'n/a';
        end
        
        if isfield(roi, 'snb')
            handles.textSignalToBackground.String = num2str(round(roi.snb,1));
        else
            handles.textSignalToBackground.String = 'n/a';
        end
        
        % should have a check for is something changed...
        if isfield(handles.data.rois, 'status')
            %if ~isfield(handles.info, 'nselected')
             %   handles.info.nselected = sum(vertcat(handles.data.rois(:,chaIndex).status) == 1);
             %   handles.info.nselected = [handles.info.nselected, handles.info.nselected];
            % end
            totalSelected = handles.info.nselected(chaIndex);
            newStr = [num2str(totalSelected), ' of ', num2str(handles.info.nrois), ' selected'];
            handles.textSelectedPercent.String = newStr;
        end
    end
end
guidata(hObject, handles);


end
