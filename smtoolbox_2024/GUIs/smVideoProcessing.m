 function varargout = smVideoProcessing(varargin)
% SMVIDEOPROCESSING MATLAB code for smVideoProcessing.fig
%      SMVIDEOPROCESSING, by itself, creates a new SMVIDEOPROCESSING or raises the existing
%      singleton*.
%
%      H = SMVIDEOPROCESSING returns the handle to a new SMVIDEOPROCESSING or the handle to
%      the existing singleton*.
%
%      SMVIDEOPROCESSING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SMVIDEOPROCESSING.M with the given input arguments.
%
%      SMVIDEOPROCESSING('Property','Value',...) creates a new SMVIDEOPROCESSING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before smVideoProcessing_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to smVideoProcessing_OpeningFcn via varargin.
%e
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help smVideoProcessing

% Last Modified by GUIDE v2.5 21-Mar-2023 11:26:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @smVideoProcessing_OpeningFcn, ...
    'gui_OutputFcn',  @smVideoProcessing_OutputFcn, ...
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


% --- Executes just before smVideoProcessing is made visible.
function smVideoProcessing_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to smVideoProcessing (see VARARGIN)

% Data shared between all GUIs
handles.data = varargin{:}; % init new fields, all here. 


% need to init all data, probs best from smHome? 


% location to load at
movegui('center')

% Choose default command line output for smVideoProcessing
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = smVideoProcessing_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
handles.data.frame{1}(1) = round(handles.slider1.Value);
guidata(hObject, handles);
guiUpdateVideo(hObject, handles, 1);


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function text_frame1_Callback(hObject, eventdata, handles)
% hObject    handle to frame1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frame1 as text
%        str2double(get(hObject,'String')) returns contents of frame1 as a double


% --- Executes during object creation, after setting all properties.
function frame1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frame1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function menu_file_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function menu_options_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function menu_plot_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function menu_navigation_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function menu_home_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function menu_traceViewer_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function menu_loadVideos_Callback(hObject, eventdata, handles)
answer = questdlg('Which Channel(s)?', 'User Answer', 'Both', 'Channel 1', 'Channel 2', 'Both');
[file, path, videos, time_s] = loadVideo();

% this is a bit wonky and should be cleaned up, especially for the case of
% loading a new video if one is already loaded

handles.data.videoPath = path;
handles.data.videoFile = file;

if ~isfield(handles.data, 'videos')
    handles.data.videos = cell(2,1);
    handles.data.videos0 = cell(2,1); % duplicate to revert back
    handles.data.time = cell(2,1);
    handles.data.frame = cell(2, 1);
end

switch answer
    case 'Both'
        handles.data.videos = videos;
        handles.data.time = time_s;
        handles.data.videos0{1} = videos{1};
        handles.data.videos0{2} = videos{2};
        
        nVideos = 2; 
    case 'Channel 1'
        handles.data.videos{1} = videos{1};
        handles.data.videos0{1} = videos{1};
        handles.data.time{1} = time_s{1};
        nVideos = 1; 
        handles.data.videos{2} =[];
        handles.data.videos0{2} = [];
        handles.data.time{2} = [];
        
    case 'Channel 2'
        handles.data.videos{2} = videos{1};
        handles.data.videos0{2} = videos{1};
        handles.data.time{2} = time_s{1};
        nVideos = 2; 
end

% gather video information (can only be 1 or two videos right now)
handles.data.nVideos = nVideos;
handles.data.videoDisplay = zeros(nVideos, 4);
handles.data.rois = [];
handles.textAOI.String = '0';

% store info
for i = 1:handles.data.nVideos
    if i == 2
        handles.data.channel2transform = 0;
    end
    if ~isempty(handles.data.videos{i})
        handles.data.frame{i} = [1;size(handles.data.videos{i},3)];
        handles.data.defaultBC{i} = [];
        % handles.data.videoDisplay(i,:) = [ handles.slider_BCmin1.Value,  handles.slider_BCmax1.Value, handles.slider_BCmin1.Value,  handles.slider_BCmax1.Value];
    end
end
guidata(hObject, handles);

% update the sliders
for i = 1:handles.data.nVideos
    if ~isempty(handles.data.videos{i})
        handles = guiInitFrameSlider(hObject, handles, i);
        handles = guiInitBCSlider(hObject, handles, i);
        guidata(hObject, handles);
    end
end

% plot the video
for i = 1:length(handles.data.videos)
    if ~isempty(handles.data.videos{i})
        % get frame, auto image BC
        guiUpdateVideo(hObject, handles, i, 1, 1);
    end
end

% --------------------------------------------------------------------
function menu_save_Callback(hObject, eventdata, handles)
guiSaveData(handles.data)

% --------------------------------------------------------------------
function menu_loadAlignment_Callback(hObject, eventdata, handles)
[file,path] = uigetfile('*.mat');
loadedtform = load([path,file]);
handles.data.alignTransform = loadedtform.alignTransform;
msgbox('Alignment transform loaded.');
guidata(hObject, handles);

% --- Executes on button press in checkbox_gauss.
function checkbox_gauss_Callback(hObject, eventdata, handles)


function edit_radius_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function edit_radius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_radius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_spacing_Callback(hObject, eventdata, handles)
% hObject    handle to edit_spacing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_spacing as text
%        str2double(get(hObject,'String')) returns contents of edit_spacing as a double


% --- Executes during object creation, after setting all properties.
function edit_spacing_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_spacing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_area_Callback(hObject, eventdata, handles)
% hObject    handle to edit_area (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_area as text
%        str2double(get(hObject,'String')) returns contents of edit_area as a double


% --- Executes during object creation, after setting all properties.
function edit_area_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_area (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in checkbox_background.
function checkbox_background_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_background (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_background

function edit_falsePositive_Callback(hObject, eventdata, handles)
% hObject    handle to edit_falsePositive (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_falsePositive as text
%        str2double(get(hObject,'String')) returns contents of edit_falsePositive as a double


% --- Executes during object creation, after setting all properties.
function edit_falsePositive_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_falsePositive (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in findAOI.
function findAOI_Callback(hObject, eventdata, handles)
% hObject    handle to findAOI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% get values
aoi.radius = round(str2double(handles.edit_radius.String));
aoi.area = round(str2double(handles.edit_area.String));
aoi.fp = round(str2double(handles.edit_falsePositive.String));
aoi.gauss = handles.checkbox_gauss.Value;
aoi.bkg = 0; 
aoi.tol = str2double(handles.editTol.String); 
aoi.method = 'GLRT';

handles.data.aoiInfo = aoi; 
guidata(hObject, handles);

% Make mean image from channel 1
if isfield(handles.data, 'videos')
    if ~isempty(handles.data.videos{1})
        s0 = 10;
        nFrames = size(handles.data.videos{1},3);
        if s0 > nFrames
            s0 = nFrames;
        end
        prompt = {'Start Frame', 'End Frame'};
        dlgtitle = 'Frame Avg.';
        dims = [1 35];
        definput = {'1', num2str(s0)};
        answer = inputdlg(prompt,dlgtitle,dims,definput);
        idx1 = str2double(answer{1});
        idx2 = str2double(answer{2});
        imageMu = mean(handles.data.videos{1}(:,:, idx1:idx2), 3);
    end
end

% channel 1
aoi.method = 'GLRT';
% handles.data.rois = findAOI(imageMu, aoi); % Otsu, GLRT
parms.method = aoi.method;
parms.radius = aoi.radius;
parms.falsePositive = aoi.fp;
parms.gaussBool = aoi.gauss;
parms.gaussTol = aoi.tol;
showResult = 0;
handles.data.rois = findAreasOfInterest(imageMu, parms); % Otsu, GLRT

handles.textAOI.String = num2str(length(handles.data.rois)); 
handles.showAOI1.Value = 1; 
guidata(hObject, handles);

% guiUpdateVideo(hObject, handles, 1);

% channel 2
if handles.data.nVideos == 2
    if isfield(handles.data, 'alignTransform')
        if isfield(handles.data.alignTransform, 'tformAOI')
            handles.data.rois = transformAOI(handles.data.rois, aoi.radius, handles.data.alignTransform);
        else
            handles.data.rois = transformAOI(handles.data.rois, aoi.radius, []);
        end
    else
        handles.data.rois = transformAOI(handles.data.rois, aoi.radius, []);
    end
    handles.showAOI2.Value = 1;
    guidata(hObject, handles);
end

for j = 1:handles.data.nVideos
    % add in centroid and bounding box for each frame
    nframe = size(handles.data.videos{j},3);
    for i = 1:length(handles.data.rois)
        handles.data.rois(i,j).Centroid = handles.data.rois(i,j).centroid;
        if isfield(handles.data.rois(i,1), 'Centroid') && ~isempty(handles.data.rois(i,j).Centroid)
            centroid = zeros(nframe,2) + handles.data.rois(i,j).Centroid;
            boundingBox = zeros(nframe,4)+ handles.data.rois(i,j).boundingBox;
            handles.data.rois(i,j).Centroid = centroid;
            handles.data.rois(i,j).boundingBox = boundingBox;
        else
            centroid = zeros(nframe,2) + handles.data.rois(i,j).centroid;
            boundingBox = zeros(nframe,4)+ handles.data.rois(i,j).boundingBox;
            handles.data.rois(i,j).Centroid = centroid;
            handles.data.rois(i,j).boundingBox = boundingBox;
        end
    end
    
    % check for out of bound (bit messy... )
    
    if min(min(handles.data.videos{j}(:,:,1))) < 1 % TEMP PATCH FOR LOADING BACKGROUND SUBTRACTED VIDEO
        idx = checkForOutofBoundROIs(handles.data.rois(:,j), handles.data.videos{j}(:,:,end)+100);
    else
        idx = checkForOutofBoundROIs(handles.data.rois(:,j), handles.data.videos{j}(:,:,end));
    end
    if ~isempty(idx)
        handles.data.rois(idx,:) = [];
    end
end

guidata(hObject, handles);
guiUpdateVideo(hObject, handles, j, 1);
if ~isempty(idx)
    handles.textAOI.String = num2str(length(handles.data.rois));
    guiUpdateVideo(hObject, handles, 1, 1);
end

% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
handles.data.frame{2}(1) = round(handles.slider2.Value);
guidata(hObject, handles);
guiUpdateVideo(hObject, handles, 2);

% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function text_frame2_Callback(hObject, eventdata, handles)
% hObject    handle to frame2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frame2 as text
%        str2double(get(hObject,'String')) returns contents of frame2 as a double


% --- Executes during object creation, after setting all properties.
function frame2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frame2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider_BCmin1_Callback(hObject, eventdata, handles)
% hObject    handle to slider_BCmin1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.slider_BCmin1.Value = round(get(hObject,'Value'));
guidata(hObject, handles);
guiUpdateVideo(hObject, handles, 1, 0, 1);


% --- Executes during object creation, after setting all properties.
function slider_BCmin1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_BCmin1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider_BCmin2_Callback(hObject, eventdata, handles)
% hObject    handle to slider_BCmin2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.slider_BCmin2.Value = round(get(hObject,'Value'));
guidata(hObject, handles);
guiUpdateVideo(hObject, handles, 2, 0, 1);

% --- Executes during object creation, after setting all properties.
function slider_BCmin2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_BCmin2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on button press in showAOI1.
function showAOI1_Callback(hObject, ~, handles)
guidata(hObject, handles);
guiUpdateVideo(hObject, handles, 1);

% --- Executes on button press in showAOI2.
function showAOI2_Callback(hObject, ~, handles)
guidata(hObject, handles);
guiUpdateVideo(hObject, handles, 2);

% --- Executes on button press in pbTimeSeries.
function pbTimeSeries_Callback(hObject, eventdata, handles)
% hObject    handle to pbTimeSeries (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isempty(handles.data.rois)
    handles.data.rois = projectAOI(handles.data.rois, handles.data.videos);
    
    % < WORK ON THIS INTEGRATION FOR SPEED INCREASE > ?
    
    % handles.data.rois = projectIntensityAOI(handles.data.rois, handles.data.videos);
    
    guidata(hObject, handles);
    guiUpdateVideo(hObject, handles, 1);
    if size(handles.data.rois,2) > 1
        guiUpdateVideo(hObject, handles, 2);
    end
    msgbox('AOIs projected');
else
    msgbox('Warning: No AOIs found to project!');
end

% Move on to image Viewer? 

% answer = questdlg('AOIs projected. Navigate to Image Viewer?', ...
% 	'Naviation', 'Yes', 'No', 'No');
% % Handle response
% switch answer
%     case 'Yes'
%         % advance to image viewer
%         
%     case 'No'
% end


function editPixels_Callback(hObject, eventdata, handles)
% hObject    handle to editPixels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editPixels as text
%        str2double(get(hObject,'String')) returns contents of editPixels as a double


% --- Executes during object creation, after setting all properties.
function editPixels_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editPixels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editTol_Callback(hObject, eventdata, handles)
% hObject    handle to editTol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editTol as text
%        str2double(get(hObject,'String')) returns contents of editTol as a double


% --- Executes during object creation, after setting all properties.
function editTol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editTol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --------------------------------------------------------------------
function menu_loadImages_Callback(hObject, eventdata, handles)
handles.data.videos = []; 
handles.data.videos.time = []; 
[path, handles.data.videos, handles.data.time] = loadMMImages();


% NEED TO UPDATE FOR SLIDER VALUE STORE -----------------------------------
% NEED TO UPDATE FOR STORING A COPY OF THE ORIGNAL VIDEO ------------------

% gather video information (can only be 1 or two videos right now)
handles.data.videoPath = path;
handles.data.videoFile = '';
handles.data.nVideos = length(handles.data.videos);
handles.data.frame = zeros(handles.data.nVideos, 2);
% handles.data.videoDisplay = zeros(handles.data.nVideos, 4);
handles.data.rois = [];
handles.textAOI.String = '0'; 

% store info
for i = 1:handles.data.nVideos
    if i == 2
        handles.data.channel2transform = 0; 
    end
    handles.data.frame(i,1) = 1;
    handles.data.frame(i,2) = size(handles.data.videos{i},3);
    image1 = handles.data.videos{i}(:,:,1);
    [mu, sigma]= normfit(image1(:));
    handles.data.videoDisplay(i,:) = [round(mu-sigma*2), round(mu+sigma*6), round(mu-sigma*2), round(mu+sigma*6)];
end
guidata(hObject, handles);

% update the slider
for i = 1:handles.data.nVideos
    handles = guiInitFrameSlider(hObject, handles, i);
    guidata(hObject, handles);
end

% plot the video
for i = 1:length(handles.data.videos)
    guiUpdateVideo(hObject, handles, i, 1);
end

% --- Executes on button press in pushbutton_drift.
function pushbutton_drift_Callback(hObject, eventdata, handles)

prompt = {'Enter Number of Frames:'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'10'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
frameIdx = round(str2double(answer{1})); 

answer = questdlg('Which Channel(s)?', 'User Answer', 'Both', 'Channel 1', 'Channel 2', 'Both');

disp('Apply drift correction...')
switch answer
    case 'Both'
        driftList = driftCorrectionFFT(handles.data.videos{1}, frameIdx);
        handles.data.rois = applyDriftFTT(handles.data.rois, driftList, [1,2]);
        guidata(hObject, handles);
        
        for j = 1:handles.data.nVideos
            guiUpdateVideo(hObject, handles, j);
            guidata(hObject, handles);
        end
        
    case 'Channel 1'
        driftList = driftCorrectionFFT(handles.data.videos{1}, frameIdx);
        handles.data.rois = applyDriftFTT(handles.data.rois, driftList, 1);
        guidata(hObject, handles);
        
    case 'Channel 2'
        driftList = driftCorrectionFFT(handles.data.videos{2}, frameIdx);
        handles.data.rois = applyDriftFTT(handles.data.rois, driftList, 2);
        guidata(hObject, handles);
        
end
msgbox('Drift Correction applied')


% --------------------------------------------------------------------
function menu_loadGlimpse_Callback(hObject, eventdata, handles)
[path, handles.data.videos, handles.data.time] = loadGlimpse();

% NEED TO UPDATE FOR SLIDE VALUE STORE ------------------------------------
% NEED TO UPDATE FOR STORING A COPY OF THE ORIGNAL VIDEO ------------------

% gather video information (can only be 1 or two videos right now)
handles.data.videoPath = path;
handles.data.videoFile = '';
handles.data.nVideos = length(handles.data.videos);
handles.data.frame = zeros(handles.data.nVideos, 2);
handles.data.videoDisplay = zeros(handles.data.nVideos, 4);
handles.data.rois = [];
handles.textAOI.String = '0'; 

% store info
for i = 1:handles.data.nVideos
    handles.data.frame(i,1) = 1;
    handles.data.frame(i,2) = size(handles.data.videos{i},3);
    image1 = handles.data.videos{i}(:,:,1);
    [mu, sigma]= normfit(image1(:));
    handles.data.videoDisplay(i,:) = [round(mu-sigma*2), round(mu+sigma*6), round(mu-sigma*2), round(mu+sigma*6)];
end
guidata(hObject, handles);

% update the slider
for i = 1:handles.data.nVideos
    handles = guiInitFrameSlider(hObject, handles, i);
    guidata(hObject, handles);
end

% plot the video
for i = 1:length(handles.data.videos)
    guiUpdateVideo(hObject, handles, i, 1);
end

% --------------------------------------------------------------------
function exportVideo_Callback(hObject, eventdata, handles)
% hObject    handle to exportVideo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
exportInfo = getExportVideoInfo(); 

bb = 0; 
scaleBar = 0; 
time_s = 0; 
umPerPixel = 0; 

if ~isempty(exportInfo)
    
    % need error check here
    n = exportInfo.channel;
    tempPath = handles.data.videoPath; 
    tempFile = ['Channel_',num2str(n),'.mp4']; % need condition for format
    tempFilePath = [tempPath, tempFile];
    
    video = handles.data.videos{n};
    time_s = handles.data.time{n}(:,end);
    nframes = size(video,3); 
    if length(time_s) ~= nframes
        prompt = {'Video Frame Rate (s)'};
        dlgtitle = 'Enter Frame Rate';
        dims = [1 35];
        definput = {'1'};
        answer = inputdlg(prompt,dlgtitle,dims,definput);
        frameRate_s = str2double(answer{1});
        time_s = (1:nframes) * frameRate_s; 
        time_s = time_s'; 
        
        %data.handles.time = {};
        %data.handles.time{1} = [1:nFrames, time_s*1000, time_s];
        %guidata(hObject, handles);
    end
        
    if exportInfo.showAOI && ~isempty(handles.data.rois)
        X = handles.data.rois(:,n);
        bb = zeros(length(X),4);
        for i=1:length(X)
            bb(i,:) = X(i).boundingBox(1,:);
        end
    else
        bb = 0;
    end
    
    if exportInfo.scaleBar_um
        scaleBar = exportInfo.scaleBar_um;
    end
    if exportInfo.umPerPixel
        umPerPixel = exportInfo.umPerPixel;
    end
    if exportInfo.frameRate_s
        frameRate_s = exportInfo.frameRate_s;
    end
    
    if n == 1
        x1 = handles.slider_BCmin1.Value;
        x2 = handles.slider_BCmax1.Value;
        
    elseif n == 2
        x1 = handles.slider_BCmin2.Value;
        x2 = handles.slider_BCmax2.Value;
    end
    
    if isnan(exportInfo.stopFrame)
        exportInfo.stopFrame = size(video,3);
    elseif  exportInfo.stopFrame >  size(video,3)
        exportInfo.stopFrame = size(video,3);
    end
    % need a few more error checks but good enough for now 7/20/23-DSW
    % write images
    writeImagesToVideo(video(:,:,exportInfo.startFrame:exportInfo.stopFrame), 'aoi', bb,...
        'time_s', time_s(exportInfo.startFrame:exportInfo.stopFrame),...
    'scaleBar', scaleBar, 'umPerPixel', umPerPixel, 'filePath', tempFilePath, ...
    'frameRate', frameRate_s, 'displayRange', [x1, x2]);
end


% --------------------------------------------------------------------
function exportFrame_Callback(hObject, eventdata, handles)
% hObject    handle to exportFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

exportInfo = getExportFrameInformation();

bb = 0; 
scaleBar = 0;
time_s = 0;
umPerPixel = 0;

if ~isempty(exportInfo)
    
    frame = exportInfo.frame;
    % need error check here
    n = exportInfo.channel;
    scaleLimits = []; 
    if n == 1
        scaleLimits = [handles.slider_BCmin1.Value, handles.slider_BCmax1.Value];
    elseif n == 2
        scaleLimits = [handles.slider_BCmin2.Value, handles.slider_BCmax2.Value];
    end
    nframes = size(handles.data.videos{n},3);
    if frame > nframes
        errordlg('Frame is out of bounds');
        return
    else
        
        tempPath = handles.data.videoPath;
        tempFile = ['Channel_',num2str(n), '-', num2str(frame)]; % need condition for format
        tempFilePath = [tempPath, tempFile];
        
        image = handles.data.videos{n}(:,:,frame);
        time_s = handles.data.time{n}(:,end);
        
        if length(time_s) ~= nframes
            prompt = {'Video Frame Rate (s)'};
            dlgtitle = 'Enter Frame Rate';
            dims = [1 35];
            definput = {'1'};
            answer = inputdlg(prompt,dlgtitle,dims,definput);
            frameRate_s = str2double(answer{1});
            time_s = frame*frameRate_s;
        else
            time_s = time_s(frame);
        end
        
        if exportInfo.showAOI && ~isempty(handles.data.rois)
            X = handles.data.rois(:,n);
            bb = zeros(length(X),4);
            for i=1:length(X)
                bb(i,:) = X(i).boundingBox(1,:);
            end
        end
        
        if exportInfo.scaleBar_um
            scaleBar = exportInfo.scaleBar_um;
        end
        
        if exportInfo.umPerPixel
            umPerPixel = exportInfo.umPerPixel;
        end
        
        % write images
        exportFrame(image, 'aoi', bb, 'time_s', time_s,...
            'scaleBar', scaleBar, 'umPerPixel', umPerPixel,...
            'scaleLimits', scaleLimits, 'filePath', tempFilePath);
    end
end


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% apply image trasnform direclty (should make into a button press at some
% point)
if handles.data.channel2transform == 0
    if isfield(handles.data, 'alignTransform')
        if isfield(handles.data.alignTransform, 'tformImage')
            % need to store refObj in align transform. Compute for now
            if ~isempty(handles.data.videos{2})
                refObj = imref2d(size(handles.data.videos{2}(:,:,1)));
                nFrames = size(handles.data.videos{2},3);
                vid2 = handles.data.videos{2};
                wbt = waitbar(0, 'Apply transform...');
                for j = 1:nFrames
                    vid2(:,:,j) = imwarp(vid2(:,:,j), refObj,...
                        handles.data.alignTransform.tformImage, 'OutputView', refObj, 'SmoothEdges', true);
                    waitbar(j/nFrames, wbt)
                end
                handles.data.videos{2} = vid2;
                handles.data.channel2transform = 1;
                guiUpdateVideo(hObject, handles, 2);
                guidata(hObject, handles);
                close(wbt);
            end
        else
            errordlg('No Alignment provided');
            
        end
    else
        errordlg('No Alignment provided')
    end
end


% --- Executes on slider movement.
function slider_BCmax1_Callback(hObject, eventdata, handles)
% hObject    handle to slider_BCmax1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.slider_BCmax1.Value = round(get(hObject,'Value'));
guidata(hObject, handles);
guiUpdateVideo(hObject, handles, 1, 0, 1);

% --- Executes during object creation, after setting all properties.
function slider_BCmax1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_BCmax1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider_BCmax2_Callback(hObject, eventdata, handles)
% hObject    handle to slider_BCmax2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.slider_BCmax2.Value = round(get(hObject,'Value'));
guidata(hObject, handles);
guiUpdateVideo(hObject, handles, 2, 0, 1);

% --- Executes during object creation, after setting all properties.
function slider_BCmax2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_BCmax2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbuttonBC1.
function pushbuttonBC1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonBC1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Auto adjust the contrast of image 1
img = handles.data.videos{1}(:,:,handles.data.frame{1}(1));
[x1,x2] = autoImageBC(img);
handles.data.defaultBC{1}(1) = x1;
handles.slider_BCmin1.Value = x1;
handles.slider_BCmax1.Value = x2;
handles.data.defaultBC{1}(2) = x2;
guidata(hObject, handles);
guiUpdateVideo(hObject, handles, 1, 0, 1);

% --- Executes on button press in pushbuttonBC2.
function pushbuttonBC2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonBC2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img = handles.data.videos{2}(:,:,handles.data.frame{2}(1));
[x1,x2] = autoImageBC(img);
handles.data.defaultBC{2}(1) = x1;
handles.slider_BCmin2.Value = x1;
handles.slider_BCmax2.Value = x2;
handles.data.defaultBC{2}(2) = x2;
guidata(hObject, handles);
guiUpdateVideo(hObject, handles, 2, 0, 1);

% --------------------------------------------------------------------
function saveVideo_Callback(hObject, eventdata, handles)
% hObject    handle to saveVideo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Option to export data back to tiff (background, drift, align)
[saveFile,savePath] = uiputfile('.tiff'); 

% hardcoded to save imagestack 2
image2 = mean(handles.data.videos{2}(:,:,100:103),3);
t = Tiff([savePath, saveFile],'w');
t.write(image2);
imwrite(image2, [savePath, saveFile])


function resetImages_Callback(hObject, eventdata, handles)
% hObject    handle to resetImages (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function resetChannel1_Callback(hObject, eventdata, handles)
% hObject    handle to resetChannel1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.data.videos{1} = handles.data.videos0{1}; 
handles = guiInitBCSlider(hObject, handles, 1, 1);
guidata(hObject, handles);
guiUpdateVideo(hObject, handles, 1, 0, 1);

% --------------------------------------------------------------------
function resetChannel2_Callback(hObject, eventdata, handles)
% hObject    handle to resetChannel2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.data.videos{2} = handles.data.videos0{2}; 
 handles.data.channel2transform = 0; 
handles = guiInitBCSlider(hObject, handles, 2, 1);
guidata(hObject, handles);
guiUpdateVideo(hObject, handles, 2, 0, 1);


% --------------------------------------------------------------------
function plotIntensityOverTime_Callback(hObject, eventdata, handles)
% hObject    handle to plotIntensityOverTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

prompt = {'Channel', 'Save Figure'};
dlgtitle = 'Enter Frame Rate';
dims = [1 35];
definput = {'1', '1'};
answer = inputdlg(prompt,dlgtitle,dims,definput);

vidIdx = str2double(answer{1});
if ~isempty(handles.data.videos{vidIdx})
    if isempty(handles.data.time{vidIdx})
        prompt = {'Video Frame Rate (s)'};
        dlgtitle = 'Enter Frame Rate';
        dims = [1 35];
        definput = {'1'};
        answer = inputdlg(prompt,dlgtitle,dims,definput);
        frameRate_s = str2double(answer{1});
        time_s = (1:nframes) * frameRate_s;
        time_s = time_s';
    else
        time_s = handles.data.time{vidIdx}(:,end); % SHOULD FIX
    end
    intensity_au = plotIntensityOverTime(handles.data.videos{vidIdx}, time_s);
else
    msg = ['Error: No Video in Channel ', answer{1}]; 
    errordlg(msg);
end

output = [handles.data.time{vidIdx}, intensity_au];
fileName = [handles.data.videoPath, 'intensityPerFrame-Channel', num2str(vidIdx),'.xlsx'];
xlswrite(fileName,output)

% --------------------------------------------------------------------
function plotAOIOverTime_Callback(hObject, eventdata, handles)
% hObject    handle to plotAOIOverTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
aoi = struct; 
aoi.radius = round(str2double(handles.edit_radius.String));
aoi.area = round(str2double(handles.edit_area.String));
aoi.fp = round(str2double(handles.edit_falsePositive.String));
aoi.gauss = handles.checkbox_gauss.Value;
aoi.gauss = 0; 
aoi.bkg = 0; 
aoi.tol = str2double(handles.editTol.String); 
aoi.method = 'GLRT';

% for now, just channel 1
channelIdx = 1; 
nFrames = size(handles.data.videos{channelIdx},3);
aoiPerFrame = zeros(nFrames,1);
wb = waitbar(0,'Find AOI in each frame...');
for i = 1:nFrames
    im = handles.data.videos{channelIdx}(:,:,i);
    aoiTemp = findAOI(im, aoi);
    aoiPerFrame(i,1) = length(aoiTemp);
    waitbar(i/nFrames,wb);
end
close(wb)

figure; 
plot(1:nFrames, aoiPerFrame,'-o')

output = [handles.data.time{channelIdx}, aoiPerFrame];
fileName = [handles.data.videoPath, 'aoiPerFrame-Channel', num2str(channelIdx),'.xlsx'];
xlswrite(fileName,output)

% --------------------------------------------------------------------
function flipImageHeader_Callback(hObject, eventdata, handles)
% hObject    handle to flipImageHeader (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function computeAndApplyDrift_Callback(hObject, eventdata, handles)
% hObject    handle to computeAndApplyDrift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function subtractBackground_Callback(hObject, eventdata, handles)
% hObject    handle to subtractBackground (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function subtractBackgroundC1_Callback(hObject, eventdata, handles)
% hObject    handle to subtractBackgroundC1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
prompt = {'Disk radius'};
dlgtitle = 'Subtract Background';
dims = [1 35];
definput = {'4'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
se = strel('disk', str2double(answer{1}));

handles.data.videos{1} = subtractBackground(handles.data.videos{1}, se);
guidata(hObject, handles);
handles = guiInitBCSlider(hObject, handles, 1, 1);
guiUpdateVideo(hObject, handles, 1, 0, 1);

% --------------------------------------------------------------------
function subtractBackgroundC2_Callback(hObject, eventdata, handles)
% hObject    handle to subtractBackgroundC2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
prompt = {'Disk radius'};
dlgtitle = 'Subtract Background';
dims = [1 35];
definput = {'4'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
se = strel('disk', str2double(answer{1}));

handles.data.videos{2} = subtractBackground(handles.data.videos{2}, se);
guidata(hObject, handles);
handles = guiInitBCSlider(hObject, handles, 2, 1);
guiUpdateVideo(hObject, handles, 2, 0, 1);

% --------------------------------------------------------------------
function flipChannel1_Callback(hObject, eventdata, handles)
% hObject    handle to flipChannel1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function flipChannel2_Callback(hObject, eventdata, handles)
% hObject    handle to flipChannel2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function flipHorizontalChannel2_Callback(hObject, eventdata, handles)
% hObject    handle to flipHorizontalChannel2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.data.videos{2} = flip(handles.data.videos{2}, 2); 
guidata(hObject, handles);
guiUpdateVideo(hObject, handles, 2, 0, 1);

% --------------------------------------------------------------------
function flipVerticalChannel2_Callback(hObject, eventdata, handles)
% hObject    handle to flipVerticalChannel2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.data.videos{2} = flip(handles.data.videos{2},1); 
guidata(hObject, handles);
guiUpdateVideo(hObject, handles, 2, 0, 1);

% --------------------------------------------------------------------
function flipHorizontalChannel1_Callback(hObject, eventdata, handles)
% hObject    handle to flipHorizontalChannel1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.data.videos{1} = flip(handles.data.videos{1},2); 
guidata(hObject, handles);
guiUpdateVideo(hObject, handles, 1, 0, 1);

% --------------------------------------------------------------------
function flipVerticalChannel1_Callback(hObject, eventdata, handles)
% hObject    handle to flipVerticalChannel1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.data.videos{1} = flip(handles.data.videos{1},1); 
guidata(hObject, handles);
guiUpdateVideo(hObject, handles, 1, 0, 1);


% --------------------------------------------------------------------
function driftCorrectionVideo_Callback(hObject, eventdata, handles)
% hObject    handle to driftCorrectionVideo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
prompt = {'Compute from channel (1 or 2): ', 'Apply to channel (0 {all}, 1, 2)', 'Frame Index (int)'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'1','0', '10'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
vidIdx = str2double(answer{1});
frameIdx = str2double(answer{3});
videoTemp = handles.data.videos{vidIdx};

switch answer{2}
    case '0'
        [newVideo, allTForm] = driftCorrectVideo(videoTemp, [], frameIdx);
        
        if vidIdx == 1
            handles.data.videos{1} = newVideo;
            if ~isempty(handles.data.videos{2})
                newVideo = driftCorrectVideo(handles.data.videos{2}, allTForm, frameIdx);
                handles.data.videos{2} = newVideo;
            end
            
        else
            handles.data.videos{2} = newVideo;
            newVideo = driftCorrectVideo(handles.data.videos{1}, allTForm, frameIdx);
            handles.data.videos{1} = newVideo;
        end
        guiUpdateVideo(hObject, handles, 1);
        guiUpdateVideo(hObject, handles, 2);
        guidata(hObject, handles);
        
    case {'1', '2'}
        [newVideo, ~] = driftCorrectVideo(videoTemp, [], frameIdx);
        handles.data.videos{vidIdx} = newVideo;
        guiUpdateVideo(hObject, handles, vidIdx);
        guidata(hObject, handles);
        
end


% --------------------------------------------------------------------
function driftCorrectionAOI_Callback(hObject, eventdata, handles)
% hObject    handle to driftCorrectionAOI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
prompt = {'Enter Number of Frames:'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'10'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
frameIdx = round(str2double(answer{1})); 

answer = questdlg('Which Channel(s)?', 'User Answer', 'Both', 'Channel 1', 'Channel 2', 'Both');

disp('Apply drift correction...')
switch answer
    case 'Both'
        driftList = driftCorrectionFFT(handles.data.videos{1}, frameIdx);
        handles.data.rois = applyDriftFTT(handles.data.rois, driftList, [1,2]);
        guidata(hObject, handles);
        
        for j = 1:handles.data.nVideos
            guiUpdateVideo(hObject, handles, j);
            guidata(hObject, handles);
        end
        
    case 'Channel 1'
        driftList = driftCorrectionFFT(handles.data.videos{1}, frameIdx);
        handles.data.rois = applyDriftFTT(handles.data.rois, driftList, 1);
        guidata(hObject, handles);
        
    case 'Channel 2'
        driftList = driftCorrectionFFT(handles.data.videos{2}, frameIdx);
        handles.data.rois = applyDriftFTT(handles.data.rois, driftList, 2);
        guidata(hObject, handles);
        
end
msgbox('Drift Correction applied')


% --------------------------------------------------------------------
function bandpassFilter_Callback(hObject, eventdata, handles)
% hObject    handle to bandpassFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function bandpassFilterC1_Callback(hObject, eventdata, handles)
% hObject    handle to bandpassFilterC1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
prompt = {'Noise Size (pixels):', 'Object Size (pixels)'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'1', '5'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
x = str2double(answer{1});
y = str2double(answer{2});
nFrames = size(handles.data.videos{1},3);
for i = 1:nFrames
    handles.data.videos{1}(:,:,i) = bpass(handles.data.videos{1}(:,:,i), x, y);
end
guidata(hObject, handles);
handles = guiInitBCSlider(hObject, handles, 1, 1);
guiUpdateVideo(hObject, handles, 1, 0, 1);

% --------------------------------------------------------------------
function bandpassFilterC2_Callback(hObject, eventdata, handles)
% hObject    handle to bandpassFilterC2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
prompt = {'Noise Size (pixels):', 'Object Size (pixels)'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'1', '5'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
x = str2double(answer{1});
y = str2double(answer{2});
nFrames = size(handles.data.videos{2},3);
for i = 1:nFrames
    handles.data.videos{2}(:,:,i) = bpass(handles.data.videos{2}(:,:,i), x, y);
end
guidata(hObject, handles);
handles = guiInitBCSlider(hObject, handles, 2, 1);
guiUpdateVideo(hObject, handles, 2, 0, 1);


% --------------------------------------------------------------------
function menuAlignment_Callback(hObject, eventdata, handles)
% hObject    handle to menuAlignment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menuAOI_Callback(hObject, eventdata, handles)
% hObject    handle to menuAOI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menuAOIParameters_Callback(hObject, eventdata, handles)
% hObject    handle to menuAOIParameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menuAlignImages_Callback(hObject, eventdata, handles)
% hObject    handle to menuAlignImages (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles.data, 'alignTransform')
    handles.data.alignTransform = struct;
    handles.data.alignTransform.tformImage = [];
    handles.data.alignTransform.alignAOI = [];
    handles.data.alignTransform.x = [];
    handles.data.alignTransform.rmse = [];
    handles.data.alignTransform.centroids = [];
end
% need error check for if video 2 event exists

image1 = handles.data.videos{1}(:,:, handles.data.frame{1}(1));
image2 = handles.data.videos{2}(:,:, handles.data.frame{2}(1));

showPlot = 1; % put in option?
[handles.data.alignTransform.tformImage, ~, refObj1, refObj2] = alignImages(image1, image2, showPlot);

% apply transform to all images in channel 2
nFrames = size(handles.data.videos{2},3);
for i = 1:nFrames
    handles.data.videos{2}(:,:,i) = imwarp(handles.data.videos{2}(:,:,i), refObj2,...
        handles.data.alignTransform.tformImage, 'OutputView', refObj1, 'SmoothEdges', true); %%%
end
guidata(hObject, handles);
handles = guiInitBCSlider(hObject, handles, 2, 1); % REALLY FIX THIS
guiUpdateVideo(hObject, handles, 2, 0, 1);


% --------------------------------------------------------------------
function menuAlignAOI_Callback(hObject, eventdata, handles)
% hObject    handle to menuAlignAOI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% grab AOI source from window (really should store in handle.data..)
if ~isfield(handles.data, 'alignTransform')
    handles.data.alignTransform = struct;
    handles.data.alignTransform.tformImage = [];
    handles.data.alignTransform.tformAOI = [];
    handles.data.alignTransform.rmse = [];
    handles.data.alignTransform.centroids = [];
end
aoi.radius = round(str2double(handles.edit_radius.String));
aoi.area = round(str2double(handles.edit_area.String));
aoi.fp = round(str2double(handles.edit_falsePositive.String));
aoi.gauss = handles.checkbox_gauss.Value;
aoi.bkg = 0; 
aoi.tol = str2double(handles.editTol.String); 
aoi.method = 'GLRT';

% ADD IN GUI FOR WHICH FRAMES TO USE

% find rois in each channel (need option for taking mean of image)
image1 = handles.data.videos{1}(:,:, handles.data.frame{1}(1));
image2 = handles.data.videos{2}(:,:, handles.data.frame{2}(1));
%image1 = mean(handles.data.videos{1}(:,:,end-3:end), 3);
%image2 = mean(handles.data.videos{2}(:,:,end-3:end), 3);

%image1 = mean(handles.data.videos{1}(:,:,100:103), 3);
%image2 = mean(handles.data.videos{2}(:,:,100:103), 3);

% get pairs
maxDist = 1;
showPlot = 1;
[tform, centroidPairs, rmse, rois1, ~] = alignAOI(image1, image2, aoi, maxDist, showPlot);
handles.data.alignTransform.tformAOI = tform;
handles.data.alignTransform.rmse = rmse;
handles.data.alignTransform.centroids = centroidPairs;

% show AOI found
handles.data.rois = rois1;
handles.data.rois = transformAOI(handles.data.rois, aoi.radius, handles.data.alignTransform);
handles.showAOI1.Value = 1;
handles.showAOI2.Value = 1;
guidata(hObject, handles);
guiUpdateVideo(hObject, handles, 1, 1, 0);
guiUpdateVideo(hObject, handles, 2, 1, 0);

% --------------------------------------------------------------------
function menuApplyImageTransform_Callback(hObject, eventdata, handles)
% hObject    handle to menuApplyImageTransform (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menyApplyAOITransform_Callback(hObject, eventdata, handles)
% hObject    handle to menyApplyAOITransform (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menuSaveAlignment_Callback(hObject, eventdata, handles)
% hObject    handle to menuSaveAlignment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles.data, 'alignTransform')
    alignTransform = handles.data.alignTransform;
    if ~isempty(alignTransform)
        path = handles.data.videoPath;
        tform_name = [path, 'alignTransform.mat'];
        [file, path] = uiputfile(tform_name);
        save([path,file], 'alignTransform');
        disp('Transform Saved:');
        disp(['  >> ', [path, file]]);
    end
end


% --- Executes on button press in pushbutton_filterAOI.
function pushbutton_filterAOI_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_filterAOI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[~, idx1] = smvideo_removeOutliers(handles.data.rois(:,1), 'maxIntensity', 0);
[~, idx2]= smvideo_removeOutliers(handles.data.rois(:,1), 'gaussSigma', 0);

% find zero pixels in last frame
img = handles.data.videos{end}(:,:,end);
if isempty(img)
    img = handles.data.videos{1}(:,:,end);
end
% check for zero values in last channel
idx3 = []; 
for i = 1:length(handles.data.rois(:,end))
    pixelList = handles.data.rois(i,end).pixelList;
    for j = 1:size(pixelList,1)
        col = pixelList(j,1);
        row = pixelList(j,2);
        if img(row, col) == 0
            idx3 = [idx3;i];
        end
    end
end

remove_index = unique([idx1; idx2; idx3]);
handles.data.rois(remove_index, :) = []; 
handles.textAOI.String = num2str(length(handles.data.rois)); 
guidata(hObject, handles);
for j = 1:2
    guiUpdateVideo(hObject, handles, j, 1);
end


% --- Executes on button press in pushbutton_coloc.
function pushbutton_coloc_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_coloc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% temporary rewrite to auto crop images in both channels from drift

% adapted from
% https://www.mathworks.com/matlabcentral/answers/297568-how-to-remove-white-background-from-image

xx = 4;
X = zeros(handles.data.nVideos,4); % xstart, xstop, ystart, ystop
if handles.data.nVideos == 1
    binaryImage = handles.data.videos{1}(:,:,end) > 0;
else
    binaryImage = handles.data.videos{2}(:,:,end) > 0;
end
binaryImageFilt = bwareafilt(binaryImage, 1);
[rows, columns] = find(binaryImageFilt);
row1 = min(rows);
row2 = max(rows);
col1 = min(columns);
col2 = max(columns);

img = binaryImageFilt(row1:row2, col1:col2);
figure; imshow(img, [])

% crop videos
for i = 1:handles.data.nVideos
    for j = 1:size(handles.data.videos{i},3)
        handles.data.videos{i}(:,:,j) =  handles.data.videos{i}(row1:row2, col1:col2, :);
    end
end
    
    

% adjust AOIs

% centroids = cell(2,1);
% for i = 1:2
%     tmp = handles.data.aoiInfo;
%     
%     parms.method = tmp.method;
%     parms.radius = tmp.radius;
%     parms.falsePositive = tmp.fp;
%     parms.falsePositive = 30;
% 
%     parms.gaussBool = tmp.gauss;
%     parms.gaussTol = tmp.tol;
%     showResult = 0;
%     
%     parms.minDist = 1;
%     nFrames = handles.data.frame{i}(end);
%     
%     if nFrames < 10
%         imageMu = mean(handles.data.videos{i}(:,:, 1:3), 3);
%     else
%         imageMu = mean(handles.data.videos{i}(:,:, 1:10), 3);
%     end
%     temp_rois = findAreasOfInterest(imageMu, parms);    % Otsu, GLRT
%     
%     % filter temp_rois
%     temp_rois = smvideo_removeOutliers(temp_rois, 'maxIntensity', 0);
%     temp_rois = smvideo_removeOutliers(temp_rois, 'gaussSigma', 0);
%     centroids{i} = vertcat(temp_rois.centroid);
% end
% 
% % Now find overlapping spots from red to green within rmse of alignment 
% maxDist = mean(handles.data.alignTransform.rmse);
% maxDist = 2;
% pairs = [];
% for i = 1:length(centroids{1})
%     dist = sqrt((centroids{1}(i,1) - centroids{2}(:,1)).^2 + (centroids{1}(i,2) - centroids{2}(:,2)).^2);
%     [mindist, c2] = min(dist);
%     if mindist <= maxDist
%         pairs = [pairs; i,c2];
%     end
% end
% 
% npairs = size(pairs,1); 
% channelNames = {'633nm', '532nm'};
% figure;
% for i = 1:2
%     img = mean(handles.data.videos{i}(:,:,1:end),3);
%     [x1,y1] = autoImageBC(img);
%     subplot(1,2,i);
%     imshow(img, 'DisplayRange', [x1,y1]);
%     hold on
%     scatter(centroids{i}(:,1), centroids{i}(:,2), 30,  'yo');
%     
%     for j = 1:npairs
%         idx = pairs(j,i);
%         scatter(centroids{i}(idx,1), centroids{i}(idx,2), 30,  'ro');
%     end
%     title(channelNames{i});
% end
% % disp(['nPairs': num2str(npairs)]);
% disp(['Total number of 633nm: ', num2str(size(centroids{1},1))]); 
% disp(['Total number of 532nm: ', num2str(size(centroids{2},1))]); 
% disp(['Total number of colocalized: ', num2str(size(pairs,1))]); 
% disp(['%% Coclocalized (ref 633): ', num2str(size(pairs,1)/ size(centroids{1},1)*100)]); 


% --------------------------------------------------------------------
function Untitled_2_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function channel1_bkg_Callback(hObject, eventdata, handles)
% hObject    handle to channel1_bkg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function channel2_bkg_Callback(hObject, eventdata, handles)
% hObject    handle to channel2_bkg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Load a background file (noise from imageJ)
[~,~,img_bkg] = loadVideo();
handles.data.videos{2} = img_bkg{1}
guidata(hObject, handles);
guiUpdateVideo(hObject, handles, 2, 0, 1);


% --------------------------------------------------------------------
function Untitled_3_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function append_channel1_Callback(hObject, eventdata, handles)
% hObject    handle to append_channel1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% load new image
[new_file, new_path, new_videos, new_time_s] = loadVideo();

% Lets assume we need to realign the new video to the first video and just
% apply it 
[tform, ~, refObj1, refObj2] = alignImages(handles.data.videos{1}(:,:,end), new_videos{1}(:,:,1), 0);
% apply tform
for i = 1:size(new_videos{1},3)
     new_videos{1}(:,:,i) = imwarp( new_videos{1}(:,:,i), refObj2,...
        tform, 'OutputView', refObj1, 'SmoothEdges', true);
end

handles.data.videos{1} = cat(3, handles.data.videos{1}, new_videos{1});
handles.data.videos0{1} = handles.data.videos{1};
handles.data.time{1} = [ handles.data.time{1}; new_time_s{1} + handles.data.time{1}(end,:)];
handles.data.frame{1} = [1; size(handles.data.time{1},1)];
guidata(hObject, handles);

% update the slider
handles = guiInitFrameSlider(hObject, handles, 1);
handles = guiInitBCSlider(hObject, handles, 1);
guidata(hObject, handles);




% --------------------------------------------------------------------
function append_channel2_Callback(hObject, eventdata, handles)
% hObject    handle to append_channel2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
