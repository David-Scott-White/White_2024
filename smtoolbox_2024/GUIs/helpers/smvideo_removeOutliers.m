function [aois, idx] = smvideo_removeOutliers(aois, filterField, showResult)
% -------------------------------------------------------------------------
% Remove the aois based on intensity values. Auto for now
% -------------------------------------------------------------------------
%
%
% Auto for now using isoutlier
% 
% David S. White 
% 2022-04-20 
% License: GPL GNU v3
% -------------------------------------------------------------------------
if nargin < 2
    showResult = 0; 
    filterField = 'sumIntensity';
end

if nargin < 3
    showResult = 0; % histogram
end

% To do: assert AOI is struct with centroid as field

if showResult
    aois0 = aois;
end
N1 = size(aois,1); 

% Main
% n1 = 1;
% n2 = 0;
% outlierMethod = 'quartile';
% remove_idx = [];
% while n1 ~= n2
%     n1 = size(aois,1);
%     switch filterField
%         case 'maxIntensity'
%             %aois(isoutlier(vertcat(aois.maxMaskIntensity),outlierMethod)) = [];
%             idx = find(isoutlier(vertcat(aois.maxMaskIntensity),outlierMethod)==1);
%         case 'avgIntensity'
%             %aois(isoutlier(vertcat(aois.avgMaskIntensity),outlierMethod)) = [];
%             idx = find(isoutlier(vertcat(aois.avgMaskIntensity),outlierMethod)==1);
%         case 'sumIntensity'
%             %aois(isoutlier(vertcat(aois.sumMaskIntensity),outlierMethod)) = [];
%             idx = find(isoutlier(vertcat(aois.sumMaskIntensity),outlierMethod)==1);
%         case 'gaussSigma'
%             %aois(isoutlier(vertcat(aois.gaussSigma),outlierMethod)) = [];  
%             idx = find(isoutlier(vertcat(aois.gaussSigma),outlierMethod)==1);
%     end
%     remove_idx = 1;
%     n2 = n1-length(idx);
%     aois(idx) = []; 
%     %n2 = size(aois, 1);
% end

% TEMP PATCH SINGLE RUN TO STORE CORRECT INDEX

outlierMethod = 'quartile';
switch filterField
    case 'maxIntensity'
        %aois(isoutlier(vertcat(aois.maxMaskIntensity),outlierMethod)) = [];
        idx = find(isoutlier(vertcat(aois.maxMaskIntensity),outlierMethod)==1);
    case 'avgIntensity'
        %aois(isoutlier(vertcat(aois.avgMaskIntensity),outlierMethod)) = [];
        idx = find(isoutlier(vertcat(aois.avgMaskIntensity),outlierMethod)==1);
    case 'sumIntensity'
        %aois(isoutlier(vertcat(aois.sumMaskIntensity),outlierMethod)) = [];
        idx = find(isoutlier(vertcat(aois.sumMaskIntensity),outlierMethod)==1);
    case 'gaussSigma'
        %aois(isoutlier(vertcat(aois.gaussSigma),outlierMethod)) = [];
        idx = find(isoutlier(vertcat(aois.gaussSigma),outlierMethod)==1);
        
end
disp(['Outliers Removed: ', num2str(length(idx)), ' | Method: ', filterField]);

if showResult
    figure;
    switch filterField
        case 'maxIntensity'
            subplot(2,1,1); histogram(vertcat(aois0.maxMaskIntensity));  title('Before');
            subplot(2,1,2); histogram(vertcat(aois.maxMaskIntensity)); title('After');
        case 'avgIntensity'
            subplot(2,1,1); histogram(vertcat(aois0.avgMaskIntensity));  title('Before');
            subplot(2,1,2); histogram(vertcat(aois.avgMaskIntensity)); title('After');
        case 'sumIntensity'
            subplot(2,1,1); histogram(vertcat(aois0.sumMaskIntensity));  title('Before');
            subplot(2,1,2); histogram(vertcat(aois.sumMaskIntensity)); title('After');
        case 'gaussSigma'
            subplot(2,1,1); histogram(vertcat(aois0.gaussSigma));  title('Before');
            subplot(2,1,2); histogram(vertcat(aois.gaussSigma)); title('After');
    end
end

