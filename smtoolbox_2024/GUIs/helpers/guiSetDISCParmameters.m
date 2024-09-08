function varargout = guiSetDISCParmameters(varargin)
% GUISETDISCPARMAMETERS MATLAB code for guiSetDISCParmameters.fig
%      GUISETDISCPARMAMETERS, by itself, creates a new GUISETDISCPARMAMETERS or raises the existing
%      singleton*.
%
%      H = GUISETDISCPARMAMETERS returns the handle to a new GUISETDISCPARMAMETERS or the handle to
%      the existing singleton*.
%
%      GUISETDISCPARMAMETERS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUISETDISCPARMAMETERS.M with the given input arguments.
%
%      GUISETDISCPARMAMETERS('Property','Value',...) creates a new GUISETDISCPARMAMETERS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before guiSetDISCParmameters_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to guiSetDISCParmameters_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help guiSetDISCParmameters

% Last Modified by GUIDE v2.5 09-Feb-2022 00:03:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @guiSetDISCParmameters_OpeningFcn, ...
                   'gui_OutputFcn',  @guiSetDISCParmameters_OutputFcn, ...
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


% --- Executes just before guiSetDISCParmameters is made visible.
function guiSetDISCParmameters_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to guiSetDISCParmameters (see VARARGIN)

% Choose default command line output for guiSetDISCParmameters
handles.output = hObject;

handles.data=struct;
handles.data.autoDISC = handles.checkboxAutoDISC.Value;
handles.data.alpha = str2double(handles.editAlpha.String);
handles.data.divisive = handles.popupmenuDivisive.String{handles.popupmenuDivisive.Value};
handles.data.agglomerative = handles.popupmenuAgglomerative.String{handles.popupmenuAgglomerative.Value};
handles.data.viterbi = handles.checkboxViterbi.Value;
handles.data.return_k = str2double(handles.editReturnKStates.String);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes guiSetDISCParmameters wait for user response (see UIRESUME)
uiwait(handles.figure1);



% --- Outputs from this function are returned to the command line.
function varargout = guiSetDISCParmameters_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% varargout{1} = handles.output;
varargout{1} = handles.data;
delete(hObject)


% --- Executes on button press in checkboxAutoDISC.
function checkboxAutoDISC_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxAutoDISC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxAutoDISC
if handles.checkboxAutoDISC.Value == 1
    handles.editAlpha.Enable = 'off';
    handles.popupmenuDivisive.Enable = 'off';
    handles.popupmenuAgglomerative.Enable = 'off';
    handles.checkboxViterbi.Enable = 'off';
    handles.editReturnKStates.Enable = 'off';
    handles.checkboxReturnKStates.Enable = 'off';
    
    handles.data.autoDISC = 1;
    handles.data.alpha = [];
    handles.data.divisive = [];
    handles.data.agglomerative = [];
    handles.data.viterbi = [];
    handles.data.return_k = [];


else
    handles.editAlpha.Enable = 'on';
    handles.popupmenuDivisive.Enable = 'on';
    handles.popupmenuAgglomerative.Enable = 'on';
    handles.checkboxViterbi.Enable = 'on';
    handles.checkboxReturnKStates.Enable = 'on';
    
    handles.data.autoDISC = 0;
    handles.data.alpha = str2double(handles.editAlpha.String);
    handles.data.divisive = handles.popupmenuDivisive.String{handles.popupmenuDivisive.Value};
    handles.data.agglomerative = handles.popupmenuAgglomerative.String{handles.popupmenuAgglomerative.Value};
    handles.data.viterbi = handles.checkboxViterbi.Value;
    handles.data.return_k = [];
    
    switch handles.checkboxReturnKStates.Value
        case 1
            handles.editReturnKStates.Enable = 'on';
            handles.data.return_k = str2double(handles.editReturnKStates.String);
        case 2
            handles.editReturnKStates.Enable = 'off';
            handles.data.return_k = 0;
    end
    
end

guidata(hObject, handles);

% --- Executes on selection change in popupmenuDivisive.
function popupmenuDivisive_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuDivisive (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuDivisive contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuDivisive
handles.data.divisive = handles.popupmenuDivisive.String{handles.popupmenuDivisive.Value};
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenuDivisive_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuDivisive (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenuAgglomerative.
function popupmenuAgglomerative_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuAgglomerative (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuAgglomerative contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuAgglomerative
handles.data.agglomerative = handles.popupmenuAgglomerative.String{handles.popupmenuAgglomerative.Value};
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenuAgglomerative_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuAgglomerative (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editAlpha_Callback(hObject, eventdata, handles)
% hObject    handle to editAlpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editAlpha as text
%        str2double(get(hObject,'String')) returns contents of editAlpha as a double
handles.data.alpha = str2double(handles.editAlpha.String);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function editAlpha_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editAlpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editReturnKStates_Callback(hObject, eventdata, handles)
% hObject    handle to editReturnKStates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editReturnKStates as text
%        str2double(get(hObject,'String')) returns contents of editReturnKStates as a double


% --- Executes during object creation, after setting all properties.
function editReturnKStates_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editReturnKStates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkboxViterbi.
function checkboxViterbi_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxViterbi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxViterbi
 handles.viterbi = handles.checkboxViterbi.Value;
guidata(hObject, handles);

% --- Executes on button press in checkboxReturnKStates.
function checkboxReturnKStates_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxReturnKStates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxReturnKStates
if get(hObject,'Value') == 1
    handles.editReturnKStates.Enable = 'on';
    handles.data.return_k = str2double(handles.editReturnKStates.String);
else
     handles.editReturnKStates.Enable = 'off';
     handles.data.return_k = 0;
end
guidata(hObject, handles);

% --- Executes on button press in pushbuttonApplySettings.
function pushbuttonApplySettings_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonApplySettings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
figure1_CloseRequestFcn(hObject, eventdata, handles)


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
% delete(hObject);
handles.data.autoDISC = handles.checkboxAutoDISC.Value;
handles.data.alpha = str2double(handles.editAlpha.String);
handles.data.divisive = handles.popupmenuDivisive.String{handles.popupmenuDivisive.Value};
handles.data.agglomerative = handles.popupmenuAgglomerative.String{handles.popupmenuAgglomerative.Value};
handles.data.viterbi = handles.checkboxViterbi.Value;
handles.data.return_k = str2double(handles.editReturnKStates.String);
guidata(hObject, handles);
uiresume()                  % resume UI which will trigger the OutputFcn
