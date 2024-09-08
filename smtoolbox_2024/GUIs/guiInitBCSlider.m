function handles = guiInitBCSlider(hObject, handles, videoIndex, recompute)

% NEED DEBUG FOR IF THE MAX VALUE IS LOWER THAN THE MIN VALUE
if nargin < 4
    recompute = 0;
end

if ~isfield(handles.data, 'defaultBC')
    handles.data.defaultBC = cell(2,1);
end
if isempty(handles.data.defaultBC{videoIndex}) || recompute
    temp = handles.data.videos{videoIndex}(:, :, 1);
    [x1, x2] = autoImageBC(temp);
    handles.data.defaultBC{videoIndex} = [x1, x2];
end

switch videoIndex
    case 1
        
        handles.slider_BCmin1.Value = handles.data.defaultBC{1}(1);
        handles.slider_BCmax1.Value = handles.data.defaultBC{1}(2);
       % handles.data.defaultBC{1} = [handles.slider_BCmin1.Value, handles.slider_BCmax1.Value];
        
        handles.slider_BCmin1.Min = 0;
        handles.slider_BCmin1.Max = handles.slider_BCmax1.Value; % IMAGE J 
        handles.slider_BCmax1.Min = handles.slider_BCmin1.Value;
        handles.slider_BCmax1.Max = 65355; % IMAGE J = 65355; more practical value for slider
        
        handles.slider_BCmin1.SliderStep = [1/(handles.slider_BCmin1.Max-1) , 1/(handles.slider_BCmin1.Max-1)];
        handles.slider_BCmax1.SliderStep = [1/(handles.slider_BCmax1.Max-1) , 1/(handles.slider_BCmax1.Max-1)];
        
        handles.text_bcmin1.String = num2str(round(handles.slider_BCmin1.Value));
        handles.text_bcmax1.String = num2str(round(handles.slider_BCmax1.Value));
        
    case 2
        handles.slider_BCmin2.Value = handles.data.defaultBC{2}(1);
        handles.slider_BCmax2.Value = handles.data.defaultBC{2}(2);
        % handles.data.defaultBC{2} = [handles.slider_BCmin2.Value, handles.slider_BCmax2.Value];
        
        handles.slider_BCmin2.Min = 0;
        handles.slider_BCmin2.Max = handles.slider_BCmax2.Value; % IMAGE J
        handles.slider_BCmax2.Min = handles.slider_BCmin2.Value;
        handles.slider_BCmax2.Max = 65355;
        
        handles.slider_BCmin2.SliderStep = [1/(handles.slider_BCmin2.Max-1) , 1/(handles.slider_BCmin2.Max-1)];
        handles.slider_BCmax2.SliderStep = [1/(handles.slider_BCmax2.Max-1) , 1/(handles.slider_BCmax2.Max-1)];
        
        handles.text_bcmin2.String = num2str(round(handles.slider_BCmin2.Value));
        handles.text_bcmax2.String = num2str(round(handles.slider_BCmax2.Value));
end