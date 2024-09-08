function handles = guiClearIdealization(hObject, handles)
chaIdx = handles.info.channelIndex;
applyToAll = handles.checkboxApplyToAll.Value;
if applyToAll
    wb = waitbar(0, ['Clearing all Idealization in Channel ', num2str(chaIdx)]);
    for i = 1:handles.info.nrois
        handles.data.rois(i, chaIdx).fit = [];
        waitbar(i/handles.info.nrois, wb);
    end
    close(wb);
else
    handles.data.rois(handles.info.roiIndex, chaIdx).fit =[];
    handles.data.rois(handles.info.roiIndex, chaIdx).events = []; 
end
guidata(hObject, handles);