function handles = guiInitFrameSlider(hObject, handles, videoIndex)
% Set default vaues for Frame Sliders

switch videoIndex
    case 1
        handles.slider1.Min = 1;
        handles.slider1.Max = handles.data.frame{videoIndex}(2);
        if handles.slider1.Max > 1
            handles.slider1.SliderStep = [1/(handles.slider1.Max-1) , 1/(handles.slider1.Max-1)];
        else
            handles.slider1.SliderStep = [0, 0];
        end
        handles.slider1.Value = 1;
        handles.text_maxFrames1.String = num2str(handles.slider1.Max);
        
    case 2
        handles.slider2.Min = 1;
        handles.slider2.Max = handles.data.frame{videoIndex}(2);
        if handles.slider2.Max > 1
            handles.slider2.SliderStep = [1/(handles.slider2.Max-1) , 1/(handles.slider2.Max-1)];
        else
            handles.slider2.SliderStep = [0, 0];
        end
        handles.slider2.Value = 1;
        handles.text_maxFrames2.String = num2str(handles.slider2.Max);
end
end
