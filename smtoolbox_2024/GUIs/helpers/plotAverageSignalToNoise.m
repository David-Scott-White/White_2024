function plotAverageSignalToNoise(rois, channel, parseBySelected)
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
SNR = vertcat(rois(:,channel).snr);

% check for existance of the figure
figName = ['Signal to Noise | Channel ', num2str(channel)];
figExist = findobj('type','figure','name',figName);
if ~isempty(figExist)
    close(figName);
end
figure('Name', figName, 'NumberTitle','off');
if parseBySelected % color selected, leave non-selected white
    selectedIndex = vertcat(rois(:,channel).status)>0; % logic
    [~, edges] = histcounts(SNR);
    histogram(SNR, edges, 'FaceColor', [0.9, 0.9, 0.9], 'EdgeColor', 'k'); hold on
    histogram(SNR(selectedIndex), edges, 'FaceColor', [0 0.4470 0.7410], 'EdgeColor', 'k');
    legend(['Unselected: n = ', num2str(length(SNR))], ['Selected: n = ', num2str(sum(selectedIndex))]);
else
    histogram(SNR, 'faceColor', [0 0.4470 0.7410]);
    legend(['n = ', num2str(length(SNR))]);
end
pbaspect([3 2 1])
ylabel('Count');
xlabel('Signal to Noise');
title(figName);

end
   