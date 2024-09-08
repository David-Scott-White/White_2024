function handles = guiIdealizationOptions(hObject, handles)
% Parse all the options for various idealization alogrithms 

if ~isfield(handles.data, 'idealParm')
    handles.data.idealParm = cell(2,1);
end

chaIdx = handles.popupmenuChannel.Value;

switch handles.popupmenuAlgorithms.String{handles.popupmenuAlgorithms.Value}
    case 'DISC'
        % have data be able to go into to save changes made
        parm = guiSetDISCParmameters();
        % quick reformat 
        parm2 = struct; 
        parm2.method = 'DISC';
        parm2.input_type = 'alpha_value';
        parm2.intput_value = parm.alpha; 
        parm2.divisive = parm.divisive;
        parm2.agglomerative = parm.agglomerative;
        parm2.viterbi = parm.viterbi; 
        parm2.return_k = parm.return_k;
        if isnan(parm2.return_k)
            parm2.return_k = 0;
        end
        handles.data.idealParm{chaIdx} = parm2; 
        
        % parse the parameters; 
        %if parm.
        %handles.data.idealParm{handles.chaIdx} = parm;
        %parm
        
        
    case 'SKM'
        
        
    case 'vbHMM'
        
    case 'Gaussian CP'
        
    case 'Threshold'
        
end
guidata(hObject, handles);