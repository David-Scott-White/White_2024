function handles = guiUpdateFilterDisplay(hObject, handles)

chaIndex = handles.info.channelIndex;
handles.editMinStates.String = num2str(handles.info.filter(chaIndex).nStates(1));
handles.editMaxStates.String = num2str(handles.info.filter(chaIndex).nStates(2));

handles.editMinSNRsig.String = num2str(handles.info.filter(chaIndex).SNRsig(1));
handles.editMaxSNRsig.String = num2str(handles.info.filter(chaIndex).SNRsig(2));

handles.editMinSNRbkg.String = num2str(handles.info.filter(chaIndex).SNRbkg(1));
handles.editMaxSNRbkg.String = num2str(handles.info.filter(chaIndex).SNRbkg(2));

if handles.info.filter(chaIndex).nStatesOn == 1
    handles.checkboxFiterStates.Value = 1;
    handles.editMinStates.Enable = 'on';
    handles.editMaxStates.Enable = 'on';
else
    handles.checkboxFiterStates.Value = 0;
    handles.editMinStates.Enable = 'off';
    handles.editMaxStates.Enable = 'off';
end

if handles.info.filter(chaIndex).SNRsigOn == 1
    handles.checkboxFilterSNRSignal.Value = 1;
    handles.editMinSNRsig.Enable = 'on';
    handles.editMaxSNRsig.Enable = 'on';
else
    handles.checkboxFilterSNRSignal.Value = 0;
    handles.editMinSNRsig.Enable = 'off';
    handles.editMaxSNRsig.Enable = 'off';
end

if handles.info.filter(chaIndex).SNRbkgOn == 1
    handles.checkboxFilterSNRbkg.Value = 1;
    handles.editMinSNRbkg.Enable = 'on';
    handles.editMaxSNRbkg.Enable = 'on';
else
    handles.checkboxFilterSNRbkg.Value = 0;
    handles.editMinSNRbkg.Enable = 'off';
    handles.editMaxSNRbkg.Enable = 'off';
end

guidata(hObject,handles); 