function rois = applyDriftFTT(rois, driftList, channels)
%% allTform and videoIndex returned by driftCorrectionFFT.m
% need to update bounding boxes
% but transform should be uniform across the images; therefore can compute
% a single drift list  (which is presumably the average of all the
% positions)
wb = waitbar(0, 'Applying drift correction...');
[numROIs, numChannels] = size(rois);
numIndex = size(driftList,1);
total = numROIs*numChannels;
iter = 0;
keepIndex = ones(numROIs, 1);
for j = 1:numChannels
    if sum(j == channels)
        for i = 1:numROIs
            roi = rois(i, j);
            centroid0 = roi.Centroid(1,:);
            bb0 = roi.boundingBox(1,:);
            diameter = bb0(4);
            for k = 1:numIndex % first index is just 0. was 2:numIndex
                centroid1 = centroid0 + driftList(k, 3:4);
                bb1  = makeBoundingBox(centroid1, diameter);
                
                % check if outside of range (will need to know image size...)
                s1 = driftList(k, 1);
                s2 = driftList(k, 2);
                rois(i,j).Centroid(s1:s2,1) = centroid1(1);
                rois(i,j).Centroid(s1:s2,2) = centroid1(2);
                rois(i,j).boundingBox(s1:s2,1)= bb1(1);
                rois(i,j).boundingBox(s1:s2,2)= bb1(2);
                rois(i,j).boundingBox(s1:s2,3)= bb1(3);
                rois(i,j).boundingBox(s1:s2,4)= bb1(4);
                
            end
            iter = iter + 1;
            waitbar(iter/total, wb);
        end
    else
        iter = numROIs+1;
        waitbar(iter/total, wb);
    end
end
close(wb);

if numROIs > 0
    figure;
    subplot(1,2,1); plot(rois(1,1).Centroid(:,1));
    subplot(1,2,2); plot(rois(1,1).Centroid(:,2));
end
