function [driftList, allTform] = driftCorrectionFFT(video, frameWidth)
% David S. White 
% 2021-09-23

% need a check for if divisible by... a bit computationally redundant rn..
numFrames = size(video,3);
if length(frameWidth)==1
    videoIndex = 1+frameWidth:frameWidth:numFrames;
    if videoIndex(end) ~= numFrames
        videoIndex = [videoIndex, numFrames];
    end
else
    videoIndex = frameWidth; 
end
nVideoIndex = length(videoIndex);
allTform = cell(nVideoIndex-1,1); 

% align each index against the first video- should work with beads
% > accounts for cummulative nature of drift over time but may be 
% less accurate as time goes on...
showPlot = 0; 
linearTransform = 0; 
image1 = video(:,:,1);
%image1 = mean(video(:,:,1:frameWidth),3);
wb = waitbar(0, 'Computing drift...');
for i = 1:length(videoIndex)-1
    image2 = video(:,:,videoIndex(i));
    %image2 = mean(video(:,:,videoIndex(i)+1:videoIndex(i+1)), 3);
    temp = alignImages(image1, image2, linearTransform, showPlot); 
    allTform{i} = temp.tform;
    waitbar(i/nVideoIndex, wb);
end
close(wb)

% or need to adjust each centroid individually after fitting to a
% similarity matrix

% Compute a dummy drift list
driftList0 = zeros(length(videoIndex), 2);
for i = 1:length(videoIndex)-1
    [x,y] = transformPointsInverse(allTform{i}, 0, 0);
    driftList0(i+1,1) = x;
    driftList0(i+1,2) = y;
end
xx = videoIndex(:);
driftList = [[1;xx(1:end-1)], [xx(1:end-1)-1 ; numFrames], driftList0(:,1), driftList0(:,2)];

% account for outliers that may have been returned
A = filloutliers(driftList(:,3),'nearest','mean');
B = filloutliers(driftList(:,4),'nearest','mean');
driftList(:,3) = A; 
driftList(:,4) = B;

% Apply a smoothing function
x = driftList(:,1); 
yy1 = smooth(x,driftList(:,3),0.5,'rloess');
yy2 = smooth(x,driftList(:,4),0.5,'rloess');
driftList(:,3) = yy1; 
driftList(:,4) = yy2;

% driftList(1,3:4) = 0; 

figure; plot(A,'-b'); hold on; plot(B,'-r'); plot(yy1,'--b', 'linewidth',1); plot(yy2, '--r', 'linewidth',1)
legend('X', 'Y', 'X Smooth', 'Y Smooth', 'location','best');
ylabel('Drift (Pixels)')
xlabel('Frame')

% cummlative application of drift list to dummy point [100,100]
% 
% % New verison. adjust for each successive image and apply cummalative
% wb = waitbar(0, 'Computing drift...');
% showPlots = 0; 
% for i = 2:length(videoIndex)
%     if i == 2
%         image1 = video(:,:,1);
%     else
%         image1 = video(:,:,i-1);
%     end
%     image2 = video(:,:,i);
%     allTform{i-1} = alignImages(image1, image2, showPlots);
%     waitbar(i/nVideoIndex, wb);
% end
% close(wb)
% 
% driftList0 = zeros(length(videoIndex), 2);
% for i = 1:length(videoIndex)-1
%     start = [100,100]; 
%     for j = 1:i
%         start = transformPointsInverse(allTform{j}, start);
%         [i, j, start];
%     end
% end


end