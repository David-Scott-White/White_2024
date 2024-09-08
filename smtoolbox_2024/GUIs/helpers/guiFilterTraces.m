function  handles = guiFilterTraces(hObject, handles)
% Filter traces based on selection criteria
% handles.info.filter

% David S. White
% 2022-02-10
% MIT
% -------------------------------------------------------------------------

% APPLY TO ALL CHANNELS OR ONLY CURRENT CHANNEL? BETTER USE FOR THE
% CHECKBOX THAN APPLY TO ALL

c = handles.info.channelIndex;
handles.info.nselected(c) = 0;
if isfield(handles.info, 'filter')
    for i = 1:handles.info.nrois
        selection = 1;
        
        if selection && handles.info.filter(c).nStatesOn
            % check number of states (ignore if not idealized)
            if handles.info.filter(c).nStatesOn
                if isfield(handles.data.rois(i,c), 'fit')
                    if ~isempty(handles.data.rois(i,c).fit)
                        nStates = size(handles.data.rois(i,c).fit.components, 1);
                        if nStates < handles.info.filter(c).nStates(1) ||  nStates > handles.info.filter(c).nStates(2)
                            selection = 0;
                        end
                    end
                end
            end
        end
        
        % check SNRsignal
        if selection && handles.info.filter(c).SNRsigOn
            snrSignal = handles.data.rois(i,c).snr;
            if snrSignal < handles.info.filter(c).SNRsig(1) || snrSignal > handles.info.filter(c).SNRsig(2)
                selection = 0;
            end
        end
        
        % checkSNR background
        if selection &&  handles.info.filter(c).SNRbkgOn
            snrBackground = handles.data.rois(i,c).snb;
            if snrBackground < handles.info.filter(c).SNRbkg(1) || snrBackground > handles.info.filter(c).SNRbkg(2)
                selection = 0;
            end
        end
        
        % store result
        
        % if apply to both channels
        handles.data.rois(i,c).status = selection;
        handles.info.nselected(c) = handles.info.nselected(c)+selection;
    end
end
