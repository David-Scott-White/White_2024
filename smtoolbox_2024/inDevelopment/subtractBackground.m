function imageBS = subtractBackground(image, se)
% Tophat filter the image using provided structural element 
% 
% requires Image Processing Toolbox

[~,~,nFrames] = size(image) ;
minN = 5;
if nFrames > minN
    wb = waitbar(0, 'Subtracting background...');
end
imageBS = image;
for i = 1:nFrames
    imageBS(:,:,i) = imtophat(image(:,:,i), se);
    if nFrames > minN
        waitbar(i/nFrames, wb);
    end
end
if nFrames > minN
    close(wb)
end
end

