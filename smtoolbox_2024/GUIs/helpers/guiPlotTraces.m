function handles = guiPlotTraces(hObject, handles, plotID, varargin)
% -------------------------------------------------------------------------
% Plot Traces in smTraceViewer
% -------------------------------------------------------------------------
%   plotID = [1,2]. Top or Bottom Plots
%
%
%
% David S. White
% 2022-02-08
% MIT
% -------------------------------------------------------------------------
%
% Notes to self, things to check for:
% which plot? Which style? Idealization or no? Option to change colors

roiIdx = handles.info.roiIndex(1);
switch plotID
    case 1
        figureHandle1 = handles.axes1;
        figureHandle2 = handles.axes3;
        plotType = handles.popupmenuView1.Value;
        
    case 2
        figureHandle1 = handles.axes2;
        figureHandle2 = handles.axes4;
        plotType = handles.popupmenuView2.Value;
        
end
cla(figureHandle1);
cla(figureHandle2);

% switch plot Type (get value here)
switch plotType
    case {1, 3}
        
        
        % Plot 1 #REORGANIZE SO CAN SWITCH PLOT2, NOT NEED TO DRAW 1 TWICE
        
        figureHandle2.Visible = 1;
        
        % get time
        if plotType == 1
            chaIdx = 1;
        else
            chaIdx = 2;
        end
        if size(handles.data.time{chaIdx},2) == 3
            time_s = handles.data.time{chaIdx}(:,3);
        else
            time_s = 1:1:length(handles.data.rois(roiIdx, chaIdx).timeSeries);
            time_s = time_s(:);
            handles.data.time{chaIdx} = [time_s, time_s,time_s];
        end
        roi = handles.data.rois(roiIdx, chaIdx);
        
        % need option for color
        cla(figureHandle1);
        plot(time_s, roi.timeSeries, '-', 'color', handles.info.channelColor(chaIdx,:), 'Parent', figureHandle1);
        set(figureHandle1, 'XLim', [0, max(time_s)])
        figureHandle1.XLabel.String =handles.info.xLabel;
        figureHandle1.YLabel.String = handles.info.ylabel;
        figureHandle1.Box = handles.info.box';
        figureHandle1.TickDir  = handles.info.TickDir;
        figureHandle1.XGrid = handles.info.grid;
        figureHandle1.YGrid = handles.info.grid;
        
        
        % Plot second plot (histogram) --> Should update for n states
        cla(figureHandle2)
        max_value = round(max(roi.timeSeries),1);
        min_value = round(min(roi.timeSeries),-1);
        if min_value == 0
            max_value = round(max(roi.timeSeries),2);
            min_value = round(min(roi.timeSeries)*-1,2)*-1;
        end
        
        if isempty(max_value) || isempty(min_value)
            return
        end
        data_range = linspace(min_value, max_value, handles.info.histogramBins); % edges
        data_counts = histcounts(roi.timeSeries, [data_range, Inf]);
        
        if length(unique(roi.timeSeries)) > 1
            data_range = linspace(min_value, max_value, handles.info.histogramBins); % edges
            data_counts = histcounts(roi.timeSeries, [data_range, Inf]);
        else
            [data_range, data_counts] = histcounts(roi.timeSeries);
        end
        
        bar(figureHandle2, data_range, data_counts,'FaceColor', handles.info.channelColor(chaIdx,:),...
            'EdgeColor', handles.info.channelColor(chaIdx,:), 'BarWidth', 1);
        hold(figureHandle2,'on');
        set(figureHandle2,'xtick',[]);
        set(figureHandle2,'ytick',[]);
        view(figureHandle2,[90,-90])
        figureHandle2.Box = 'off';
        
        % set the axis to match the time series plot (alignaxes) axis.
        set(figureHandle2, 'XLim',get(figureHandle1,'Ylim'));
        
        % controls for if ideall...
        
        if isfield(roi, 'fit') && ~isempty(roi.fit)
            
            % Is Highlight an Event?  MAKE INFO handles.info
            if handles.checkboxAdjustEvents.Value && handles.info.channelIndex == chaIdx
                
                cla(figureHandle1);
                event = roi.events(handles.info.currentEvent,:);
                s1 = event(1);
                s2 = event(2);
                s0 = roi.fit.components(event(4),2);
                plot(time_s, roi.timeSeries, '-', 'color', handles.info.backgroundColor, 'Parent', figureHandle1);
                hold(figureHandle1,'on');
                plot(time_s, roi.fit.ideal, '-', 'color', handles.info.backgroundIdealColor,...
                    'linewidth', 1, 'Parent', figureHandle1);
                
                plot(time_s(s1:s2), roi.timeSeries(s1:s2),'-', 'linewidth', 1, 'color', handles.info.channelColor(chaIdx,:),'Parent', figureHandle1);
                plot(time_s(s1:s2), zeros(event(3),1)+s0, '-', 'linewidth', handles.info.highhightWidth,...
                    'color', handles.info.highlightColor,'Parent', figureHandle1);
                
                hold(figureHandle1,'off');
                set(figureHandle1, 'XLim', [0, max(time_s)])
                figureHandle1.XLabel.String =handles.info.xLabel;
                figureHandle1.YLabel.String = handles.info.ylabel;
                figureHandle1.Box = handles.info.box';
                figureHandle1.TickDir  = handles.info.TickDir;
                figureHandle1.XGrid = handles.info.grid;
                figureHandle1.YGrid = handles.info.grid;
                
            else
                
                hold(figureHandle1,'on');
                
                % add lines for each state
                for i = 1:size(roi.fit.components,1)
                    yline(roi.fit.components(i,2),'--', 'color', [0.7,0.7,0.7],  'Parent', figureHandle1);
                end

                plot(time_s, roi.fit.ideal, '-', 'color', handles.info.idealColor,...
                    'linewidth', handles.info.idealWidth, 'Parent', figureHandle1);
                
                hold(figureHandle1,'off');
                
                
            end
            
            % Gaussian Plot
            components = roi.fit.components;
            n_components = size(components,1);
            
            % Evalute and plot each component individually
            hold(figureHandle1,'on');
            gauss_fit_all = zeros(size(data_range));
            for n = 1:n_components
                w = components(n,1);       % weight
                mu = components(n,2);      % mu
                sigma = components(n,3);   % sigma
                
                % Evaluate each gaussian distribtution
                norm_dist_pdf = normpdf(data_range, mu, sigma).*trapz(data_range, data_counts);
                
                % store sum for gauss_fit_all
                gauss_fit_all = gauss_fit_all + w .* normpdf(data_range, mu, sigma);
                
                % convert PDF to distribution
                norm_dist = norm_dist_pdf * round(w,2);
                
                % plot this gaussian component onto the histogram
                plot(figureHandle2, data_range, norm_dist,  '--', 'color',  handles.info.idealColor, 'linewidth', handles.info.idealWidth-0.5);
                
            end
            % Compute the sum of all Gaussians
            gauss_fit_all = gauss_fit_all.* trapz(data_range, data_counts);
            plot(figureHandle2, data_range, gauss_fit_all, '-', 'color', handles.info.idealColor, 'linewidth', handles.info.idealWidth);
            hold(figureHandle1,'off');
        end
        
    case {2, 4}
        
        
        % plot AOI
        cla(figureHandle1);
        cla(figureHandle2);
        % get time
        if plotType == 2
            chaIdx = 1;
            
        else
            chaIdx = 2;
        end
        roi = handles.data.rois(roiIdx,chaIdx);
        figureHandle2.Visible = 0;
        
        if isfield(roi, 'spot')
            % temp patchhhhhhh *******************************************
            if ~isfield(roi, 'events')  || isempty(roi.events)
                spots = roi.spot;
                % get width of figure handle
                %figureWidth = floor(figureHandle1.Position(3));
                %figureHeight = floor(figureHandle1.Position(4));
                %[spotWidth, spotHeight, nSpots] = size(roi.spot);
                %nSpotsWidth = floor(figureWidth/spotWidth);
                %nSpotsHeight =floor(figureHeight/spotHeight);
                
                % need a way to control the size..
                %width = handles.info.plotSpotWidth;
                
                %if width > length(spots)
                %    width = length(spots);
                %end
                %width = 100;
                %nSpotsTemp = nSpots;
                %while rem(nSpotsTemp,width)
                %    nSpotsTemp = nSpotsTemp+1;
                %end
                %[mu, sigma] = normfit(spots(:));
                %mag = 1;
                % [nSpotsTemp/width width];
                %out = imtile(spots, 'GridSize',[nSpotsTemp/width width], 'thumbnailSize', [spotWidth*mag, spotHeight*mag]);
                % out = imtile(spots, 'GridSize',[5,2], 'thumbnailSize', [spotWidth*mag, spotHeight*mag]);
                % out = imtile(spots, 'thumbnailSize', [spotWidth*mag, spotHeight*mag]);
                
                
                im = roi.spot;
                [mu,sigma] = normfit(im(:));
                
                if size(im,3) > 10
                    %nImages = size(im,3);
                    nImages = ceil(size(im,3)/10)*10;
                    n1 = 1:ceil(sqrt(nImages));
                    n2 = ceil(nImages./n1);
                    n3 = n2./n1;
                    figureHandle1.Units = 'pixels';
                    w = figureHandle1.Position(3);
                    h = figureHandle1.Position(4);
                    figureHandle1.Units = 'normalized';
                    aspectratio = w/h;
                    [~,i] = min(abs(aspectratio-n3));
                    ncol = n2(i);
                    nrow = n1(i);
                else
                    ncol = size(im,3);
                    nrow = 1;
                end
                im2 = imtile(im, 'GridSize', [nrow, ncol]);

                if mu == 0
                    imshow(im2,  'DisplayRange', [0,1], 'Parent', figureHandle1);
                else
                    imshow(im2,  'DisplayRange', [mu-sigma, mu+6*sigma], 'Parent', figureHandle1);
                end
                
                % need to figure out ideal
                % need to figure out how to make it fill the entire plot...
                if isfield(roi, 'fit') && ~isempty(roi.fit)
                    
                    hold(figureHandle1,'on');
                    
                    arrayMap = reshape(1:nSpotsTemp, [width, round(nSpotsTemp)/width ]);
                    arrayMap = arrayMap';
                    ss = handles.data.rois(roiIdx, chaIdx).fit.ideal;
                    
                    % need case for if all idealized...
                    boundIdx = find(ss>min(ss));
                    nX = zeros(length(boundIdx),2);
                    
                    for i = 1:length(boundIdx)
                        j = boundIdx(i);
                        [jy, jx] = find(arrayMap==j);
                        nX(i,1) = jx;
                        nX(i,2) = jy;
                    end
                    for i = 1:length(boundIdx)
                        xPos = nX(i,1)*spotWidth*mag-spotWidth*mag + 0.5;
                        yPos = nX(i,2)*spotHeight*mag-spotHeight*mag + 0.5;
                        rectangle('Parent', figureHandle1, 'Position',...
                            [xPos, yPos, (spotWidth)*mag, (spotHeight)*mag], 'EdgeColor',...
                            handles.info.idealPlotSpotColor, 'lineWidth', handles.info.idealPlotSpotWidth);
                    end
                    
                    hold(figureHandle1,'off');
                end
                
            else
                events = roi.events; 
                nEvents = size(events,1);
                [spotWidth, spotHeight, ~] = size(roi.spot);
                
                
                imageMu = zeros(spotWidth, spotHeight, nEvents); 
                for k = 1:nEvents
                    s1 = events(k,1); 
                    s2 = events(k,2);
                    imageMu(:,:,k) = mean(roi.spot(:,:,s1:s2),3);
                end
                [mu, sigma] = normfit(imageMu(:));
                
                
                if nEvents > 5
                    %nImages = size(im,3);
                    %nImages = ceil(size(nEvents,3)/10)*10;
                    n1 = 1:ceil(sqrt(nEvents));
                    n2 = ceil(nEvents./n1);
                    n3 = n2./n1;
                    cla(figureHandle1);
                    figureHandle1.Units = 'pixels';
                    w = figureHandle1.Position(3);
                    h = figureHandle1.Position(4);
                    figureHandle1.Units = 'normalized';
                    aspectratio = w/h;
                    [~,i] = min(abs(aspectratio-n3));
                    ncol = n2(i);
                    nrow = n1(i);
                else
                    ncol = size(imageMu,3);
                    nrow = 1;
                end
                out = imtile(imageMu, 'GridSize',[nrow ,ncol], 'ThumbnailSize', [50,50]);
                if mu == 0
                    imshow(out,  'DisplayRange', [0,1], 'Parent', figureHandle1, 'InitialMagnification', 10);
                else
                    imshow(out,  'DisplayRange', [mu-sigma, mu+5*sigma], 'Parent', figureHandle1,'InitialMagnification', 100);
                end
            end
            
        end
        
end
guidata(hObject, handles);


end



