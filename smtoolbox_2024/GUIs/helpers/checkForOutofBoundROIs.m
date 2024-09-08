function zeroIdx = checkForOutofBoundROIs(rois, image)
% David S. White 
% 2022-02-03
% MIT
zeroIdx = [];
[X, ~, ~] = size(image);
pixelList = find(image==0); 
if ~isempty(pixelList)
    zeroIdx = []; 
    for i = 1:length(rois)
        px = rois(i).pixelList;
        for j = 1:length(px)
            p = X*(px(j,1)-1) + px(j,2);
            if sum(p == pixelList)
                zeroIdx = [zeroIdx; i]; 
               break
            end
        end
    end
end