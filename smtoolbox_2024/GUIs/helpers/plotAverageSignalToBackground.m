function plotAverageSignalToBackground(rois, channel, parseBySelected)
%% Plot average Signal to background (ie snr_bkg) of all the traces 
%
% function called from smTraceViewer.m
%
% input = rois (from data.rois). 
%
% David S. White 
% 2022-02-11
% MIT 
% -------------------------------------------------------------------------
SNB = vertcat(rois(:,channel).snb);

% check for existance of the figure
figName = ['Signal to Background | Channel ', num2str(channel)];
figExist = findobj('type','figure','name',figName);
if ~isempty(figExist)
    close(figName);
end
figure('Name', figName, 'NumberTitle','off');
if parseBySelected % color selected, leave non-selected white
    selectedIndex = vertcat(rois(:,channel).status)>0; % logic
    [~, edges] = histcounts(SNB);
    counts1 = histcounts(SNB, edges); 
    counts2 =  histcounts(SNB(selectedIndex), edges); 
    histogram(SNB, edges, 'FaceColor', [0.9, 0.9, 0.9], 'EdgeColor', 'k'); hold on
    histogram(SNB(selectedIndex), edges, 'FaceColor', [0 0.4470 0.7410], 'EdgeColor', 'k');
    legend(['Unselected: n = ', num2str(length(SNB))], ['Selected: n = ', num2str(sum(selectedIndex))]);
else
    histogram(SNB, 'faceColor', [0 0.4470 0.7410]);
    legend(['n = ', num2str(length(SNB))]);
end
pbaspect([3 2 1])
ylabel('Count');
xlabel('Signal to Background');
title(figName);

end
   