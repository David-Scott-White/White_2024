function varargout = smHome(varargin)
% SMHOME MATLAB code for smHome.fig
%      SMHOME, by itself, creates a new SMHOME or raises the existing
%      singleton*.
%
%      H = SMHOME returns the handle to a new SMHOME or the handle to
%      the existing singleton*.
%
%      SMHOME('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SMHOME.M with the given input arguments.
%
%      SMHOME('Property','Value',...) creates a new SMHOME or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before smHome_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to smHome_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help smHome

% Last Modified by GUIDE v2.5 02-Sep-2021 14:17:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @smHome_OpeningFcn, ...
                   'gui_OutputFcn',  @smHome_OutputFcn, ...
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


% --- Executes just before smHome is made visible.
function smHome_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to smHome (see VARARGIN)

% Data shared between all GUIs
handles.data = varargin{:};
movegui('center')

% Choose default command line output for smHome
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes smHome wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = smHome_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_save.
function pushbutton_save_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton_export.
function pushbutton_export_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_export (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton_append.
function pushbutton_append_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_append (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton_load.
function pushbutton_load_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton_align.
function pushbutton_align_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_align (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(findobj('Name','alignChannels'));
close(findobj('Name','smHome'));
alignChannels(handles.data);

% --- Executes on button press in pushbutton_image.
function pushbutton_image_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(findobj('Name','smVideoProcessing'));
close(findobj('Name','smHome'));
smVideoProcessing(handles.data);

% --- Executes on button press in pushbutton_trace.
function pushbutton_trace_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_trace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton_analysis.
function pushbutton_analysis_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_analysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton_simulation.
function pushbutton_simulation_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_simulation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
