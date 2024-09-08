function handles = guiAdjustNumStates(hObject, handles, p)
% Manually increase or decrease the number of states
% p = -1 (decrease state) or 1 (increase state)

roiIndex = handles.info.roiIndex;
channelIndex = handles.info.channelIndex;
if isfield(handles.data.rois(roiIndex, channelIndex), 'fit')
    if ~isempty(handles.data.rois(roiIndex, channelIndex).fit)
        nstates = size(handles.data.rois(roiIndex, channelIndex).fit.components,1);
        minstates = 1;
        maxstates = size(handles.data.rois(roiIndex,channelIndex).fit.all_ideal,2);
        k = nstates+p;
        if k >= minstates && k <= maxstates
            ss = handles.data.rois(roiIndex, channelIndex).fit.all_ideal(:,k);
            ts = handles.data.rois(roiIndex, channelIndex).timeSeries;
             [viterbi_data_fit] = runViterbi(ts, ss,3);
            [comps, ideal, class] = computeCenters(handles.data.rois(roiIndex, channelIndex).timeSeries,...
                viterbi_data_fit);
            handles.data.rois(roiIndex, channelIndex).fit.components = comps;
            handles.data.rois(roiIndex, channelIndex).fit.ideal = ideal;
            handles.data.rois(roiIndex, channelIndex).fit.class = class;
            handles.data.rois(roiIndex, channelIndex).fit.parameters.return_k = k;
            handles.data.rois(roiIndex, channelIndex).events = ...
                findEvents(handles.data.rois(roiIndex, channelIndex).fit.class);
        end
        guidata(hObject, handles);
    end
end