function handles = guiIdealizeTrace(hObject, handles)
% Idealize the trace(s)
% see guiIdealizationOptions.m

% chaIdx = handles.popupmenuChannel.Value;
chaIdx = handles.info.channelIndex;
if ~isfield(handles.data, 'idealParm') || isempty(handles.data.idealParm{chaIdx})
        % set default
        handles = guiIdealizationOptions(hObject, handles);
        guidata(hObject, handles);
end

parm = handles.data.idealParm{chaIdx};
if parm.viterbi == 1
    parm.viterbi = 3;
end
applyToAll = handles.checkboxApplyToAll.Value;

switch handles.data.idealParm{chaIdx}.method
    
    case 'DISC'
        if applyToAll
            % modify waitbar for some other options (cancel, show nrois)
            wb = waitbar(0, ['DISC: Idealizing all | Channel ', num2str(chaIdx)]);
            for i = 1:handles.info.nrois
                handles.data.rois(i, chaIdx).fit = ...
                    runDISC(handles.data.rois(i, chaIdx).timeSeries, parm);
                
                % store events
                handles.data.rois(i,chaIdx).events =...
                    findEvents(handles.data.rois(i, chaIdx).fit.class); 
                
                waitbar(i/handles.info.nrois, wb);
            end
            close(wb);
        else
            roiIdx = handles.info.roiIndex;
            handles.data.rois(roiIdx, chaIdx).fit = ...
                runDISC(handles.data.rois(roiIdx, chaIdx).timeSeries, parm);
            
            % store events
            handles.data.rois(roiIdx,chaIdx).events =...
                    findEvents(handles.data.rois(roiIdx, chaIdx).fit.class); 
        end
              
    case 'SKM'
        states = localGetStatesSKM();
        if applyToAll
            wb = waitbar(0, ['vbFRET: Idealizing all | Channel ', num2str(chaIdx)]);
            for i = 1:handles.info.nrois
                [components,ideal, class] = runSKM(handles.data.rois(i, chaIdx).timeSeries, states);
                handles.data.rois(i,chaIdx).fit.components = components;
                handles.data.rois(i,chaIdx).fit.ideal = ideal;
                handles.data.rois(i,chaIdx).fit.class = class;
                handles.data.rois(i,chaIdx).fit.metrics = [];
                handles.data.rois(i,chaIdx).fit.all_ideal = [];
                handles.data.rois(i,chaIdx).fit.parameters = 'SKM';
                
                % store events
                handles.data.rois(i,chaIdx).events =...
                    findEvents(handles.data.rois(i, chaIdx).fit.class);
                
                waitbar(i/handles.info.nrois, wb);
            end
            close(wb);
        else
            roiIdx = handles.info.roiIndex;
            [components,ideal, class] = runSKM(handles.data.rois(roiIdx, chaIdx).timeSeries, states);
            handles.data.rois(roiIdx,chaIdx).fit.components = components;
            handles.data.rois(roiIdx,chaIdx).fit.ideal = ideal;
            handles.data.rois(roiIdx,chaIdx).fit.class = class;
            handles.data.rois(roiIdx,chaIdx).fit.metrics = [];
            handles.data.rois(roiIdx,chaIdx).fit.all_ideal = [];
            handles.data.rois(roiIdx,chaIdx).fit.parameters = 'SKM';
        end

        
    case 'vbHMM'
        if applyToAll
            % modify waitbar for some other options (cancel, show nrois)
            wb = waitbar(0, ['vbFRET: Idealizing all | Channel ', num2str(chaIdx)]);
            for i = 1:handles.info.nrois
                [components,ideal, class,outF,allIdeal] = runVBFRET(handles.data.rois(i, chaIdx).timeSeries);
                handles.data.rois(i,chaIdx).fit.components = components;
                handles.data.rois(i,chaIdx).fit.ideal = ideal;
                handles.data.rois(i,chaIdx).fit.class = class;
                handles.data.rois(i,chaIdx).fit.metrics = outF;
                handles.data.rois(i,chaIdx).fit.all_ideal = allIdeal;
                handles.data.rois(i,chaIdx).fit.parameters = [];
                
                % store events
                handles.data.rois(i,chaIdx).events =...
                    findEvents(handles.data.rois(i, chaIdx).fit.class);
                
                waitbar(i/handles.info.nrois, wb);
            end
            close(wb);
        else
            roiIdx = handles.info.roiIndex;
            [components,ideal, class,outF,allIdeal] = runVBFRET(handles.data.rois(roiIdx, chaIdx).timeSeries);
            handles.data.rois(roiIdx,chaIdx).fit.components = components;
            handles.data.rois(roiIdx,chaIdx).fit.ideal = ideal;
            handles.data.rois(roiIdx,chaIdx).fit.class = class;
            handles.data.rois(roiIdx,chaIdx).fit.metrics = outF;
            handles.data.rois(roiIdx,chaIdx).fit.all_ideal = allIdeal;
            handles.data.rois(roiIdx,chaIdx).fit.parameters = [];
            
            % store events
            handles.data.rois(roiIdx,chaIdx).events =...
                findEvents(handles.data.rois(roiIdx, chaIdx).fit.class);
        end
        
        
    case 'Gaussian CP'
        
    case 'Threshold'
        roi = handles.data.rois(handles.info.roiIndex, handles.info.channelIndex);
        ts = roi.timeSeries;
        
        prompt = {'Bin Frames above:', 'Viterbi (boolean):'};
        dlgtitle = 'Threshold Input';
        dims = [1 35];
        definput = {num2str(round(mean(ts)*1.1)),'1'};
        answer = inputdlg(prompt,dlgtitle,dims,definput);
        
        ss = ones(length(ts),1);
        ss(ts >= str2double(answer(1))) = 2;
        c0 = computeCenters(ts, ss);
        % refine within std of bin
        x = 2;
        ss(ss >= (c0(2,2)-c0(2,3)*x)) = 2;
        [comps, ideal, class] = computeCenters(ts, ss);
        
        if str2double(answer(2))
            data_fit = runViterbi(ts, ss, 10);
            [comps, ideal, class] = computeCenters(ts, data_fit);
        end
        handles.data.rois(handles.info.roiIndex, handles.info.channelIndex).fit.components = comps;
        handles.data.rois(handles.info.roiIndex, handles.info.channelIndex).fit.ideal = ideal;
        handles.data.rois(handles.info.roiIndex, handles.info.channelIndex).fit.class = class;
        handles.data.rois(handles.info.roiIndex, handles.info.channelIndex).fit.method = 'threshold';
        handles.data.rois(handles.info.roiIndex, handles.info.channelIndex).events =...
                    findEvents(handles.data.rois(handles.info.roiIndex, handles.info.channelIndex).fit.class); 
        guidata(hObject, handles);
        
    case 'Single CP' % CHANGE to PELT in future
        
        prompt = {'Max Change Points:', 'Penalty:'};
        dlgtitle = 'PELT Input';
        dims = [1 35];
        definput = {'1','BIC_GMM'};
        answer = inputdlg(prompt,dlgtitle,dims,definput);
        
        maxChangePoints = str2double(answer{1});
        clusterMethod = answer{2}; 
        
        if applyToAll
              wb = waitbar(0, ['Single CP: Idealizing all | Channel ', num2str(chaIdx)]);
            for i = 1:handles.info.nrois
                roi = handles.data.rois(i, handles.info.channelIndex);
                ipt = findchangepts(roi.timeSeries, 'MaxNumChanges',maxChangePoints,'Statistic', 'mean');
                ipt =[1; ipt; length(roi.timeSeries)];
                ss = roi.timeSeries*0;
                x = length(ipt)-1;
                
                for k = 1:length(ipt)-1
                    s1 = ipt(k);
                    s2 = ipt(k+1);
                    ss(s1:s2) = x;
                    x = x-1;
                end
                
                [~, ideal] = computeCenters(roi.timeSeries, ss);
                all_ideal = aggCluster(roi.timeSeries, ideal);
                [metrics, best_fit] = computeIC(roi.timeSeries, all_ideal,clusterMethod, 1);
                data_fit = all_ideal(:,best_fit);
                
                [comps, ideal, class] = computeCenters(roi.timeSeries, data_fit);
                handles.data.rois(i, handles.info.channelIndex).fit.components = comps;
                handles.data.rois(i, handles.info.channelIndex).fit.ideal = ideal;
                handles.data.rois(i, handles.info.channelIndex).fit.class = class;
                handles.data.rois(i, handles.info.channelIndex).fit.all_ideal = all_ideal;
                handles.data.rois(i, handles.info.channelIndex).fit.method = 'threshold';
                handles.data.rois(i, handles.info.channelIndex).events =...
                    findEvents(handles.data.rois(i, handles.info.channelIndex).fit.class);
                waitbar(i/handles.info.nrois, wb);
            end
            close(wb)
        else
            
            roi = handles.data.rois(handles.info.roiIndex, handles.info.channelIndex);
            ipt = findchangepts(roi.timeSeries, 'MaxNumChanges',maxChangePoints,'Statistic', 'mean');
            
            ipt =[1; ipt; length(roi.timeSeries)];
                
            ss = roi.timeSeries*0;
            x = length(ipt)-1;
            
            for k = 1:length(ipt)-1
                s1 = ipt(k);
                s2 = ipt(k+1); 
                ss(s1:s2) = x; 
                x = x-1; 
            end
            [~, ideal] = computeCenters(roi.timeSeries, ss);
            all_ideal = aggCluster(roi.timeSeries, ideal);
            [metrics, best_fit] = computeIC(roi.timeSeries, all_ideal,clusterMethod, 1);
            data_fit = all_ideal(:,best_fit);
                
            [comps, ideal, class] = computeCenters(roi.timeSeries, data_fit);
            handles.data.rois(handles.info.roiIndex, handles.info.channelIndex).fit.components = comps;
            handles.data.rois(handles.info.roiIndex, handles.info.channelIndex).fit.ideal = ideal;
            handles.data.rois(handles.info.roiIndex, handles.info.channelIndex).fit.class = class;
            handles.data.rois(handles.info.roiIndex, handles.info.channelIndex).fit.all_ideal = all_ideal;
            handles.data.rois(handles.info.roiIndex, handles.info.channelIndex).fit.method = 'threshold';
            handles.data.rois(handles.info.roiIndex, handles.info.channelIndex).events =...
                findEvents(handles.data.rois(handles.info.roiIndex, handles.info.channelIndex).fit.class);
        end
        guidata(hObject, handles);
        
end
guidata(hObject, handles);


% local functinos 
    function states = localGetStatesSKM()
        prompt = {'K States (int) or State Guesses (comma):'};
        dlgtitle = 'Threshold Input';
        dims = [1 45];
        definput = {'2'};
        answer = inputdlg(prompt,dlgtitle,dims,definput);
        answer = answer{1};
        if contains(answer,',')
            % states guesses provided
            % remove spaces
            answer = strrep(answer, ' ', '');
            % split by comma
            states = str2double(strsplit(answer,','));
            states = states(:);
        else
            states = str2double(answer);
        end
    end
end