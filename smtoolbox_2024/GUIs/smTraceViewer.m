function varargout = smTraceViewer(varargin)
% SMTRACEVIEWER MATLAB code for smTraceViewer.fig
%      SMTRACEVIEWER, by itself, creates a new SMTRACEVIEWER or raises the existing
%      singleton*.
%
%      H = SMTRACEVIEWER returns the handle to a new SMTRACEVIEWER or the handle to
%      the existing singleton*.
%
%      SMTRACEVIEWER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SMTRACEVIEWER.M with the given input arguments.
%
%      SMTRACEVIEWER('Property','Value',...) creates a new SMTRACEVIEWER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before smTraceViewer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to smTraceViewer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help smTraceViewer

% Last Modified by GUIDE v2.5 10-Feb-2022 21:55:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @smTraceViewer_OpeningFcn, ...
    'gui_OutputFcn',  @smTraceViewer_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before smTraceViewer is made visible.
function smTraceViewer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to smTraceViewer (see VARARGIN)

% Choose default command line output for smTraceViewer
handles.output = hObject;

% Data shared between all GUIs
% handles.data = varargin{:}; % init new fields, all here.
if ~isempty(varargin)
    handles.data = varargin{:};
    handles = guiInitSMTraceViewer(hObject, handles);
    guidata(hObject, handles);
end
handles = guiSetInitialParamtersTraceViewer(hObject, handles); % handles.info
handles = guiUpdateFilterDisplay(hObject, handles); % set defaults in window

% location to load at
movegui('center')

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes smTraceViewer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = smTraceViewer_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2

% --------------------------------------------------------------------
function menu_file_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_load_Callback(hObject, eventdata, handles)
% hObject    handle to menu_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file,path] = uigetfile('.mat');
if ~file
    return
end

temp = load([path,file]);
if isfield(temp, 'data')
    handles.data = temp.data;
    handles.data.aoiPath = path;
    handles.data.aoiFile = file;
    handles = guiInitSMTraceViewer(hObject, handles);
    guidata(hObject, handles);
else
    % have some error
end

% --------------------------------------------------------------------
function menu_save_Callback(hObject, eventdata, handles)
% hObject    handle to menu_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles, 'data')
    if isfield(handles.data, 'rois')
        temppath = [handles.data.aoiPath, handles.data.aoiFile];
        [savefile, savepath] = uiputfile(temppath);
        data = handles.data;
        disp('Saving...')
        save([savepath,savefile], 'data', '-v7.3');
        disp(['  >> Saved:', [savepath, savefile]]);
    end
end


% --- Executes on button press in pushbuttonSetOptions.
function pushbuttonSetOptions_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSetOptions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guiIdealizationOptions(hObject, handles);
guidata(hObject, handles);

% store in the correct channel

% --- Executes on selection change in popupmenuAlgorithms.
function popupmenuAlgorithms_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuAlgorithms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuAlgorithms contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuAlgorithms
chaIdx = handles.info.channelIndex;
handles.data.idealParm{chaIdx}.method = handles.popupmenuAlgorithms.String{handles.popupmenuAlgorithms.Value};
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenuAlgorithms_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuAlgorithms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkboxApplyToAll.
function checkboxApplyToAll_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxApplyToAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxApplyToAll


% --- Executes on button press in pushbuttonIdealize.
function pushbuttonIdealize_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonIdealize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guiIdealizeTrace(hObject, handles);

guidata(hObject, handles);

guiPlotTraces(hObject, handles, 1); % should be smarter to only redraw one
guiPlotTraces(hObject, handles, 2);
guiUpdateTraceInfo(hObject, handles);

% --------------------------------------------------------------------
function menuPlot_Callback(hObject, eventdata, handles)
% hObject    handle to menuPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menuStateOccupancy_Callback(hObject, eventdata, handles)
% hObject    handle to menuStateOccupancy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menuPlotDwellTimes_Callback(hObject, eventdata, handles)
% hObject    handle to menuPlotDwellTimes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
roi = handles.data.rois(handles.info.roiIndex, :);
figure;
plot(roi(1).timeSeries); hold on; plot(roi(2).timeSeries);

% --------------------------------------------------------------------
function menuPlotSNR_Callback(hObject, eventdata, handles)
% hObject    handle to menuPlotSNR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.info.nselected(handles.info.channelIndex) ~= handles.info.nrois
    parseBySelected = 1;
else
    parseBySelected = 0;
end
plotAverageSignalToNoise(handles.data.rois, handles.info.channelIndex, parseBySelected);

% --------------------------------------------------------------------
function menuPlotSBR_Callback(hObject, eventdata, handles)
% hObject    handle to menuPlotSBR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.info.nselected(handles.info.channelIndex) ~= handles.info.nrois
    parseBySelected = 1;
else
    parseBySelected = 0;
end
plotAverageSignalToBackground(handles.data.rois, handles.info.channelIndex, parseBySelected);

% --- Executes on selection change in popupmenu4.
function popupmenu4_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu4


% --- Executes during object creation, after setting all properties.
function popupmenu4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenuView2.
function popupmenuView2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuView2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuView2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuView2
guiPlotTraces(hObject, handles, 2);

% --- Executes during object creation, after setting all properties.
function popupmenuView2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuView2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenuChannel.
function popupmenuChannel_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuChannel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuChannel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuChannel
handles.info.channelIndex = handles.popupmenuChannel.Value;
guidata(hObject, handles);
handles = guiUpdateTraceInfo(hObject, handles);
guidata(hObject, handles);
handles = guiUpdateFilterDisplay(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenuChannel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuChannel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenuView1.
function popupmenuView1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuView1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuView1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuView1
guiPlotTraces(hObject, handles, 1);

% --- Executes during object creation, after setting all properties.
function popupmenuView1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuView1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on slider movement.
function sliderAOI_Callback(hObject, eventdata, handles)
% hObject    handle to sliderAOI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.info.currentEvent = 1;
handles.info.roiIndex(1) = round(handles.sliderAOI.Value);
handles.sliderAOI.Value = round(handles.sliderAOI.Value);
handles.textCurrentAOI.String = num2str(handles.sliderAOI.Value);
guiPlotTraces(hObject, handles, 1);
guiPlotTraces(hObject, handles, 2);
guiUpdateTraceInfo(hObject, handles);

% --- Executes during object creation, after setting all properties.
function sliderAOI_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderAOI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on key press with focus on figure1 or any of its controls.
function figure1_WindowKeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

switch eventdata.Key
    case 'leftarrow'
        handles.info.roiIndex = handles.info.roiIndex-1;
        handles.info.currentEvent = 1;
        if handles.info.roiIndex <= 0
            handles.info.roiIndex = 1;
        end
        handles.sliderAOI.Value =  handles.info.roiIndex;
        handles.textCurrentAOI.String = num2str(handles.sliderAOI.Value);
        guiPlotTraces(hObject, handles, 1);
        guiPlotTraces(hObject, handles, 2);
        guiUpdateTraceInfo(hObject, handles);
        
    case 'rightarrow'
        handles.info.roiIndex = handles.info.roiIndex+1;
        handles.info.currentEvent = 1;
        if handles.info.roiIndex > handles.info.nrois
            handles.info.roiIndex = handles.info.nrois;
        end
        handles.sliderAOI.Value =  handles.info.roiIndex;
        handles.textCurrentAOI.String = num2str(handles.sliderAOI.Value);
        guiPlotTraces(hObject, handles, 1);
        guiPlotTraces(hObject, handles, 2);
        guiUpdateTraceInfo(hObject, handles);
        
    case 'uparrow'
        if  handles.data.rois(handles.info.roiIndex, handles.info.channelIndex).status == 0
            for i = 1:handles.info.nchannels
                handles.data.rois(handles.info.roiIndex, i).status = 1;
            end
            handles.info.nselected = handles.info.nselected+1;
        end
        % guiUpdateTraceInfo(hObject, handles); % should replace this with the specific update.
        handles.textStatus.String = '0';
        %totalSelected = handles.info.nselected(handles.info.channelIndex);
        %newStr = [num2str(totalSelected), ' of ', num2str(handles.info.nrois), ' selected'];
        %handles.textSelectedPercent.String = newStr;
        guidata(hObject, handles);
        
        
        %         if ~isfield(handles.data.rois, 'status') % need init function
        %             [handles.data.rois.status] = deal(0); % default is selected
        %         end
        %         roiIdx = handles.info.roiIndex(1);
        %         if  handles.data.rois(roiIdx, handles.info.channelIndex).status == 0
        %             for i = 1:handles.info.nchannels
        %                 handles.data.rois(roiIdx, i).status = 1; % assign to all channels
        %             end
        %             handles.info.nselected =   [handles.info.nselected + 1, handles.info.nselected + 1];
        %
        %         end
        % guiUpdateTraceInfo(hObject, handles);
        handles.textStatus.String = '1';
        %totalSelected = handles.info.nselected(handles.info.channelIndex);
        % handles.textSelectedPercent.String = [num2str(totalSelected), ' of ', num2str(handles.info.nrois), ' selected'];
        %handles.textSelectedPercent.String = num2str(totalSelected);
        %guidata(hObject, handles);
        
    case  'downarrow' % NEED BUTTON TO APPLY TO ALL CHANNELS OR ONLY CURRENT
        if  handles.data.rois(handles.info.roiIndex, handles.info.channelIndex).status == 1
            for i = 1:handles.info.nchannels
                handles.data.rois(handles.info.roiIndex, i).status = 0;
            end
            handles.info.nselected = handles.info.nselected-1;
        end
        % guiUpdateTraceInfo(hObject, handles); % should replace this with the specific update.
        handles.textStatus.String = '0';
        %totalSelected = handles.info.nselected(handles.info.channelIndex);
        %newStr = [num2str(totalSelected), ' of ', num2str(handles.info.nrois), ' selected'];
        %handles.textSelectedPercent.String = newStr;
        guidata(hObject, handles);
        
    case 'space'
        j = handles.popupmenuChannel.Value +    1;
        if j > handles.info.nchannels
            j = 1;
        end
        handles.popupmenuChannel.Value = j;
        handles.info.channelIndex = j;
        guiUpdateTraceInfo(hObject, handles);
        
    case {'a', 'shift', 'm'}
        handles = guiIdealizeTrace(hObject, handles);
        guidata(hObject, handles);
        guiPlotTraces(hObject, handles, 1); % should be smarter to only redraw one
        guiPlotTraces(hObject, handles, 2);
        guiUpdateTraceInfo(hObject, handles);
        
    case {'c'}
        handles = guiClearIdealization(hObject, handles);
        guidata(hObject, handles);
        guiPlotTraces(hObject, handles, 1); % should be smarter to only redraw one
        guiPlotTraces(hObject, handles, 2);
        guiUpdateTraceInfo(hObject, handles);
        
    case 'period'
        % find next selected traces (status > 0)
        if isfield(handles, 'data') % would make sense to store a list of all? update as selected...
            i = handles.info.roiIndex;
            j = handles.info.channelIndex;
            done = 0;
            while ~done
                i = i + 1;
                if i <= size(handles.data.rois,1)
                    if handles.data.rois(i,j).status > 0
                        handles.info.roiIndex = i;
                        handles.info.currentEvent = 1;
                        handles.sliderAOI.Value =  handles.info.roiIndex;
                        handles.textCurrentAOI.String = num2str(handles.sliderAOI.Value);
                        guiPlotTraces(hObject, handles, 1);
                        guiPlotTraces(hObject, handles, 2);
                        guiUpdateTraceInfo(hObject, handles);
                        guidata(hObject, handles);
                        done = 1;
                    end
                else
                    done = 1;
                    % no updates
                end
            end
        end
        
    case 'comma'
        if isfield(handles, 'data') % would make sense to store a list of all? update as selected...
            i = handles.info.roiIndex;
            j = handles.info.channelIndex;
            done = 0;
            while ~done
                i = i - 1;
                if i > 0
                    if handles.data.rois(i,j).status > 0
                        handles.info.roiIndex = i;
                        handles.info.currentEvent = 1;
                        handles.sliderAOI.Value =  handles.info.roiIndex;
                        handles.textCurrentAOI.String = num2str(handles.sliderAOI.Value);
                        guiPlotTraces(hObject, handles, 1);
                        guiPlotTraces(hObject, handles, 2);
                        guiUpdateTraceInfo(hObject, handles);
                        guidata(hObject, handles);
                        done = 1;
                    end
                else
                    % no more selected, no updates
                    done = 1;
                end
            end
        end
        
                    
    case {'w','x'}
        %find next unselected traces (status == 0)
        if isfield(handles, 'data') % would make sense to store a list of all? update as selected...
            i = handles.info.roiIndex;
            j = handles.info.channelIndex;
            done = 0;
            while ~done
                i = i + 1;
                if i <= size(handles.data.rois,1)
                    if handles.data.rois(i,j).status == 0
                        handles.info.roiIndex = i;
                        handles.info.currentEvent = 1;
                        handles.sliderAOI.Value =  handles.info.roiIndex;
                        handles.textCurrentAOI.String = num2str(handles.sliderAOI.Value);
                        guiPlotTraces(hObject, handles, 1);
                        guiPlotTraces(hObject, handles, 2);
                        guiUpdateTraceInfo(hObject, handles);
                        guidata(hObject, handles);
                        done = 1;
                    end
                else
                    done = 1;
                    % no updates
                end
            end
        end
        
    case {'q','z'}
        % find next unselected traces (status == 0)
        if isfield(handles, 'data') % would make sense to store a list of all? update as selected...
            i = handles.info.roiIndex;
            j = handles.info.channelIndex;
            done = 0;
            while ~done
                i = i - 1;
                if i <= size(handles.data.rois,1)
                    if handles.data.rois(i,j).status == 0
                        handles.info.roiIndex = i;
                        handles.info.currentEvent = 1;
                        handles.sliderAOI.Value =  handles.info.roiIndex;
                        handles.textCurrentAOI.String = num2str(handles.sliderAOI.Value);
                        guiPlotTraces(hObject, handles, 1);
                        guiPlotTraces(hObject, handles, 2);
                        guiUpdateTraceInfo(hObject, handles);
                        guidata(hObject, handles);
                        done = 1;
                    end
                else
                    done = 1;
                    % no updates
                end
            end
        end
        
end


% --- Executes on button press in pushbuttonStatesUp.
function pushbuttonStatesUp_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonStatesUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guiAdjustNumStates(hObject, handles, 1);
guidata(hObject,handles);
guiPlotTraces(hObject, handles, 1); % should be smarter to only redraw one
guiPlotTraces(hObject, handles, 2);
guiUpdateTraceInfo(hObject, handles);

% --- Executes on button press in pushbuttonStatesDown.
function pushbuttonStatesDown_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonStatesDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guiAdjustNumStates(hObject, handles, -1);
guidata(hObject,handles);
guiPlotTraces(hObject, handles, 1); % should be smarter to only redraw one
guiPlotTraces(hObject, handles, 2);
guiUpdateTraceInfo(hObject, handles);

% --- Executes on button press in pushbuttonToggleForward.
function pushbuttonToggleForward_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonToggleForward (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.checkboxAdjustEvents.Value == 1
    i = handles.info.roiIndex;
    j = handles.info.channelIndex;
    if isfield(handles.data.rois(i,j), 'events')
        k =  handles.info.currentEvent + 1;
        if k <= size(handles.data.rois(i,j).events,1)
            handles.info.currentEvent = k;
            guidata(hObject, handles);
            guiPlotTraces(hObject, handles, 1);
            guiPlotTraces(hObject, handles, 2);
            guiUpdateTraceInfo(hObject, handles);
        end
    end
end

% --- Executes on button press in pushbuttonToggleReverse.
function pushbuttonToggleReverse_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonToggleReverse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.checkboxAdjustEvents.Value == 1
    i = handles.info.roiIndex;
    j = handles.info.channelIndex;
    if isfield(handles.data.rois(i,j), 'events')
        k =  handles.info.currentEvent - 1;
        if k > 0
            handles.info.currentEvent = k;
            % which plots to update (could make into a function, do both for now...)
            guidata(hObject, handles);
            guiPlotTraces(hObject, handles, 1);
            guiPlotTraces(hObject, handles, 2);
            guiUpdateTraceInfo(hObject, handles);
        end
    end
end

% --- Executes on button press in pushbuttonEventStateUp.
function pushbuttonEventStateUp_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonEventStateUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guiManualAdjustStates(hObject, handles, 1);
guidata(hObject, handles);

% --- Executes on button press in pushbuttonEventStateDown.
function pushbuttonEventStateDown_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonEventStateDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guiManualAdjustStates(hObject, handles, -1);
guidata(hObject, handles);

% --- Executes on button press in checkboxAdjustEvents.
function checkboxAdjustEvents_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxAdjustEvents (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxAdjustEvents
if handles.checkboxAdjustEvents.Value == 1
    if isfield(handles, 'data')
        if isfield(handles.data, 'rois')
            if isfield(handles.data.rois, 'fit')
                if ~isempty(handles.data.rois(handles.info.roiIndex, handles.info.channelIndex).fit)
                    handles.textToggleEvents.Enable = 'on';
                    handles.pushbuttonToggleReverse.Enable = 'on';
                    handles.pushbuttonToggleForward.Enable = 'on';
                    handles.textToggleEventState.Enable = 'on';
                    handles.pushbuttonEventStateDown.Enable = 'on';
                    handles.pushbuttonEventStateUp.Enable = 'on';
                    handles.textModifyNumStates.Enable = 'off';
                    handles.pushbuttonStatesDown.Enable = 'off';
                    handles.pushbuttonStatesUp.Enable = 'off';
                    
                    % prevent possible errors
                    guidata(hObject,handles);
                    guiPlotTraces(hObject, handles, 1); % should be smarter to only redraw one
                    guiPlotTraces(hObject, handles, 2);
                    guiUpdateTraceInfo(hObject, handles);
                else
                    handles.checkboxAdjustEvents.Value = 0;
                end
                
                
            else
                handles.checkboxAdjustEvents.Value = 0;
            end
        else
            handles.checkboxAdjustEvents.Value = 0;
        end
    else
        handles.checkboxAdjustEvents.Value = 0;
    end
else
    handles.textToggleEvents.Enable = 'off';
    handles.pushbuttonToggleReverse.Enable = 'off';
    handles.pushbuttonToggleForward.Enable = 'off';
    handles.textToggleEventState.Enable = 'off';
    handles.pushbuttonEventStateDown.Enable = 'off';
    handles.pushbuttonEventStateUp.Enable = 'off';
    handles.textModifyNumStates.Enable = 'on';
    handles.pushbuttonStatesDown.Enable = 'on';
    handles.pushbuttonStatesUp.Enable = 'on';
    guidata(hObject,handles);
    guiPlotTraces(hObject, handles, 1); % should be smarter to only redraw one
    guiPlotTraces(hObject, handles, 2);
    guiUpdateTraceInfo(hObject, handles);
end


% --- Executes on button press in pushbuttonApplyFilter.
function pushbuttonApplyFilter_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonApplyFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guiFilterTraces(hObject, handles);
guidata(hObject, handles);
guiUpdateTraceInfo(hObject, handles);

% --- Executes on button press in pushbuttonIntersectChannels.
function pushbuttonIntersectChannels_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonIntersectChannels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of pushbuttonIntersectChannels
if handles.info.nchannels > 1
    for i = 1:handles.info.nrois
        if sum(horzcat(handles.data.rois(i,:).status)) ~= handles.info.nchannels
            [handles.data.rois(i,:).status] = deal(0);
        end
        x = sum(vertcat(handles.data.rois(:,1).status)==1);
        handles.info.nselected = [x,x];
    end
end
guidata(hObject, handles);
guiUpdateTraceInfo(hObject, handles);

% --- Executes on button press in pushbuttonClearIdeal.
function pushbuttonClearIdeal_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonClearIdeal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guiClearIdealization(hObject, handles);
guidata(hObject, handles);
guiPlotTraces(hObject, handles, 1); % should be smarter to only redraw one
guiPlotTraces(hObject, handles, 2);
guiUpdateTraceInfo(hObject, handles);

% --------------------------------------------------------------------
function menuNavigation_Callback(hObject, eventdata, handles)
% hObject    handle to menuNavigation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menuAdvanced_Callback(hObject, eventdata, handles)
% hObject    handle to menuAdvanced (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menuNavigationHome_Callback(hObject, eventdata, handles)
% hObject    handle to menuNavigationHome (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menuNavigationImageProcessing_Callback(hObject, eventdata, handles)
% hObject    handle to menuNavigationImageProcessing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menuNavigationAnalysis_Callback(hObject, eventdata, handles)
% hObject    handle to menuNavigationAnalysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menuPlotAvgIntensity_Callback(hObject, eventdata, handles)
% hObject    handle to menuPlotAvgIntensity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.info.nselected(handles.info.channelIndex) ~= handles.info.nrois
    parseBySelected = 1;
else
    parseBySelected = 0;
end
plotAverageIntensity(handles.data.rois, handles.info.channelIndex, parseBySelected);



function editMinStates_Callback(hObject, eventdata, handles)
% hObject    handle to editMinStates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editMinStates as text
%        str2double(get(hObject,'String')) returns contents of editMinStates as a double
i = handles.info.channelIndex;
handles.info.filter(i).nStates(1) = str2double(handles.editMinStates.String);
guidata(hObject,handles);
handles = guiUpdateFilterDisplay(hObject, handles);

% --- Executes during object creation, after setting all properties.
function editMinStates_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMinStates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkboxFiterStates.
function checkboxFiterStates_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxFiterStates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxFiterStates
i = handles.info.channelIndex;
if handles.checkboxFiterStates.Value == 1
    handles.editMinStates.Enable = 'on';
    handles.editMaxStates.Enable = 'on';
    handles.info.filter(i).nStatesOn = 1;
    handles.info.filter(i).nStates(1) = str2double(handles.editMinStates.String);
    handles.info.filter(i).nStates(2) = str2double(handles.editMaxStates.String);
else
    handles.info.filter(i).nStatesOn = 0;
    handles.editMinStates.Enable = 'off';
    handles.editMaxStates.Enable = 'off';
end
guidata(hObject,handles);
handles = guiUpdateFilterDisplay(hObject, handles);

function editMaxStates_Callback(hObject, eventdata, handles)
% hObject    handle to editMaxStates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editMaxStates as text
%        str2double(get(hObject,'String')) returns contents of editMaxStates as a double
i = handles.info.channelIndex;
handles.info.filter(i).nStates(2) = str2double(handles.editMaxStates.String);
guidata(hObject,handles);
handles = guiUpdateFilterDisplay(hObject, handles);

% --- Executes during object creation, after setting all properties.
function editMaxStates_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMaxStates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editMinSNRsig_Callback(hObject, eventdata, handles)
% hObject    handle to editMinSNRsig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editMinSNRsig as text
%        str2double(get(hObject,'String')) returns contents of editMinSNRsig as a double
i = handles.info.channelIndex;
handles.info.filter(i).SNRsig(1) = str2double(handles.editMinSNRsig.String);
guidata(hObject,handles);
handles = guiUpdateFilterDisplay(hObject, handles);

% --- Executes during object creation, after setting all properties.
function editMinSNRsig_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMinSNRsig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkboxFilterSNRSignal.
function checkboxFilterSNRSignal_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxFilterSNRSignal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxFilterSNRSignal
i = handles.info.channelIndex;
if handles.checkboxFilterSNRSignal.Value == 1
    handles.editMinSNRsig.Enable = 'on';
    handles.editMaxSNRsig.Enable = 'on';
    handles.info.filter(i).SNRsigOn = 1;
    handles.info.filter(i).SNRsig(1) = str2double(handles.editMinSNRsig.String);
    handles.info.filter(i).SNRsig(2) = str2double(handles.editMaxSNRsig.String);
    
else
    handles.info.filter(1).SNRsigOn = 0;
    handles.editMinSNRsig.Enable = 'off';
    handles.editMaxSNRsig.Enable = 'off';
    guidata(hObject,handles);
end
guidata(hObject,handles);
handles = guiUpdateFilterDisplay(hObject, handles);

function editMaxSNRsig_Callback(hObject, eventdata, handles)
% hObject    handle to editMaxSNRsig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editMaxSNRsig as text
%        str2double(get(hObject,'String')) returns contents of editMaxSNRsig as a double
i = handles.info.channelIndex;
handles.info.filter(i).SNRsig(2) = str2double(handles.editMaxSNRsig.String);
guidata(hObject,handles);
handles = guiUpdateFilterDisplay(hObject, handles);

% --- Executes during object creation, after setting all properties.
function editMaxSNRsig_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMaxSNRsig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editMinSNRbkg_Callback(hObject, eventdata, handles)
% hObject    handle to editMinSNRbkg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editMinSNRbkg as text
%        str2double(get(hObject,'String')) returns contents of editMinSNRbkg as a double
i = handles.info.channelIndex;
handles.info.filter(i).SNRbkg(1) = str2double(handles.editMinSNRbkg.String);
guidata(hObject,handles);
handles = guiUpdateFilterDisplay(hObject, handles);

% --- Executes during object creation, after setting all properties.
function editMinSNRbkg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMinSNRbkg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkboxFilterSNRbkg.
function checkboxFilterSNRbkg_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxFilterSNRbkg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxFilterSNRbkg
i = handles.info.channelIndex;
if handles.checkboxFilterSNRbkg.Value == 1
    handles.editMinSNRbkg.Enable = 'on';
    handles.editMaxSNRbkg.Enable = 'on';
    handles.info.filter(i).SNRbkg(1) = str2double(handles.editMinSNRbkg.String);
    handles.info.filter(i).SNRbkg(2) = str2double(handles.editMaxSNRbkg.String);
    handles.info.filter(i).SNRbkgOn = 1;
else
    handles.info.filter(i).SNRbkgOn = 0;
    handles.editMinSNRbkg.Enable = 'off';
    handles.editMaxSNRbkg.Enable = 'off';
end
guidata(hObject,handles);
handles = guiUpdateFilterDisplay(hObject, handles);

function editMaxSNRbkg_Callback(hObject, eventdata, handles)
% hObject    handle to editMaxSNRbkg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editMaxSNRbkg as text
%        str2double(get(hObject,'String')) returns contents of editMaxSNRbkg as a double
i = handles.info.channelIndex;
handles.info.filter(i).SNRbkg(2) = str2double(handles.editMaxSNRbkg.String);
guidata(hObject,handles);
handles = guiUpdateFilterDisplay(hObject, handles);

% --- Executes during object creation, after setting all properties.
function editMaxSNRbkg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMaxSNRbkg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonSelectAll.
function pushbuttonSelectAll_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSelectAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
for j = 1:handles.info.nchannels
    for i = 1:handles.info.nrois
        handles.data.rois(i,j).status = 1;
    end
end
handles.info.nselected = [handles.info.nrois, handles.info.nrois];
guidata(hObject,handles);
handles = guiUpdateTraceInfo(hObject, handles);
guidata(hObject,handles);



function textCurrentAOI_Callback(hObject, eventdata, handles)
% hObject    handle to textCurrentAOI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of textCurrentAOI as text
%        str2double(get(hObject,'String')) returns contents of textCurrentAOI as a double
i = round(str2double(handles.textCurrentAOI.String));
handles.info.roiIndex = i;
handles.sliderAOI.Value =  handles.info.roiIndex;
handles.textCurrentAOI.String = num2str(handles.sliderAOI.Value);
guiPlotTraces(hObject, handles, 1);
guiPlotTraces(hObject, handles, 2);
guiUpdateTraceInfo(hObject, handles);


% --------------------------------------------------------------------
function menuAdvDispSetting_Callback(hObject, eventdata, handles)
% hObject    handle to menuAdvDispSetting (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
