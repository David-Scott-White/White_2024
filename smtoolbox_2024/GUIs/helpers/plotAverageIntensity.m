function plotAverageIntensity(rois, channel, parseBySelected)
%% Plot average intensity of all the traces 

% function called from smTraceViewer.m
% 
% Need option for if selected vs not selected 
%
% input = rois (from data.rois). 
%
% David S. White 
% 2022-02-09
% MIT 
% -------------------------------------------------------------------------
[nrois, ~] = size(rois);
avgIntensity = zeros(nrois,1);
for i = 1:nrois
    avgIntensity(i,1) = mean(rois(i,channel).timeSeries(:,end));
end
% check for existance of the figure
figName = ['Average Intensity | Channel ', num2str(channel)];
figExist = findobj('type','figure','name',figName);
if ~isempty(figExist)
    close(figName);
end
figure('Name', figName, 'NumberTitle','off');
if parseBySelected % color selected, leave non-selected white
    selectedIndex = vertcat(rois(:,channel).status)>0; % logic
    % plot all 
    [~, edges] = histcounts(avgIntensity);
    histogram(avgIntensity, edges, 'FaceColor', [0.9, 0.9, 0.9], 'EdgeColor', 'k'); hold on
    histogram(avgIntensity(selectedIndex), edges, 'FaceColor', [0 0.4470 0.7410], 'EdgeColor', 'k');
    legend(['Unselected: n = ', num2str(length(avgIntensity))], ['Selected: n = ', num2str(sum(selectedIndex))]);
    
else
    histogram(avgIntensity, 'faceColor', [0 0.4470 0.7410]);
    legend(['n = ', num2str(length(avgIntensity))]);
end
% color by selected vs not selected
pbaspect([3 2 1])
ylabel('Count')
xlabel('Average Intensity')
set(gca,'yscale','log')
title(figName);

end
   