function y = plotIntensityOverTime(video, time_s, scale)

if ~exist('scale', 'var') || isempty(scale)
    scale = 0; 
end

[~,~, nFrames] = size(video); 
y = zeros(nFrames,1); 
for i = 1:nFrames
    im = video(:,:,i); 
    im = imcrop(im, [100, 100, 300, 300]);
    y(i) = mean(im(:));
end

hh = figure('Name', 'Intensity Over Time'); 
if ~scale
    plot(time_s, y);
    xlabel('Time (s)');
    ylabel('Intensity (au)');
else    
    y = (y-min(y))./(max(y)-min(y));
    plot(time_s, y);
    xlabel('Time (s)');
    ylabel('Relative Intensity');
    
end


