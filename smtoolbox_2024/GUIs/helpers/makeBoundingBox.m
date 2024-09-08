function [boundingBox, pixelList] = makeBoundingBox(centroids, diameter)

if diameter == 2
    radius = 1; 
else
    radius = (diameter-1)/2;
    if radius > 1 && mod(radius,2)
        radius = radius + 0.5;
    end
end

numCentroids = size(centroids,1);
boundingBox = zeros(numCentroids, 4);
boundingBox(:,1) = floor(centroids(:,1))-radius; 
boundingBox(:,2) = floor(centroids(:,2))-radius; 
boundingBox(:,3:4) = radius*2+1;

end