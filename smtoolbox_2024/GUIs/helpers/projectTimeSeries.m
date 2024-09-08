function rois = projectTimeSeries(rois, videos, time_s, nPixels, subtractBackground)
% David S. White
% 2021-09-19

% Will need to edit for drift correction....

% Arguments
if nargin < 5 || isempty(nPixels)
    nPixels = 4;
end

[nrois, nchannels] = size(rois);

% Projection
if subtractBackground
    R = 15;
    H = 5;
    se=strel('ball', R, H);
end

% code in a wait bar, option for background subtraction
keep = ones(nrois,1);
for j = 1:nchannels
    nFrames = size(videos{j},3);
    h = waitbar(0, ['Projecting Channel ', num2str(j)]);
    for k = 1:nFrames
        if subtractBackground
            image = videos{j}(:,:,k)-imopen( videos{j}(:,:,k), se);
        else
            image = videos{j}(:,:,k);
        end
        
        for i = 1:nrois
            
           % if i==87 && j == 2
           %     disp('Stop')
           % end

            % Note, cropping by bounding box adds and extra dimension for some reason..
            if k == 1
                rois(i,j).spotInt = zeros(rois(i,j).boundingBox(k,3)+1, rois(i,j).boundingBox(k,4)+1, nFrames);
                rois(i,j).spot = zeros(7,7,nFrames);
                rois(i,j).timeSeries = zeros(nFrames,3); % frame time intensity
                rois(i,j).pixels = zeros(9,2,nFrames);
            end
            
            %%
            % ideally, store this info elsewhere to save memory. fine for now
            rois(i,j).timeSeries(:,1) = time_s{j}(:,1);
            rois(i,j).timeSeries(:,2) = time_s{j}(:,end);
            bb = rois(i,j).boundingBox(k,:);
            bb(3:4) = bb(3:4);
            spotInt = imcrop(image, bb);
            
            % figure; imshow(spotInt, [],'InitialMagnification','fit')
            % pixelList = boundBoxToPixels(bb);
            
            
            %%
            bb2  = bb; 
            bb2(1) = bb2(1)-1; 
            bb2(2) = bb2(2)-1; 
            bb2(3:4) = 6;
            spot = imcrop(image,bb2); 
            %figure; imshow(spot, [],'InitialMagnification','fit')
            
            if size(spot) == size(rois(i,j).spot, 1:2)
                rois(i,j).spot(:,:,k) = spot;
                 rois(i,j).spotInt(:,:,k) = spotInt;
                % rois(i,j).timeSeries(k,3) = sum(sort(spot(:),'descend'),1:nPixels);
                % rois(i,j).timeSeries(k,3) = sum(spotInt(:));
                % spotInt =spotInt;
                
                %spotInt(:,1) = []; 
                %spotInt(:, end) = [];
                %spotInt(1,:) = []; 
                %spotInt(end,:) = [];
              
                 rois(i,j).timeSeries(k,3) = sum(spotInt(:));
%                intSum = 0;
%                for mm = 1:size(pixelList,1)
%                    intSum = intSum + image(pixelList(mm,2), pixelList(mm,1));
%                end
%                rois(i,j).timeSeries(k,3) = intSum; 
            else
                keep(i) = 0;
            end
            % rois(i,j).timeSeries(:,3) = rois(i,j).timeSeries(:,3)-min(rois(i,j).timeSeries(:,3));
        end
        waitbar(k/nFrames, h);
    end
    close(h);
end
rois = rois(keep==1,:);

end