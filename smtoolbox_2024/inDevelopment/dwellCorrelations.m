function [h, r, p] = dwellCorrelations(x, y, binWidth, showPoints, minB, maxB)
% -------------------------------------------------------------------------
% Dwell Time Correaltion 
% -------------------------------------------------------------------------
% Summary: 
% 
% Input:
%   dwellPairs = [N,2] matrix of x-values and y-values
% 
% Output:
%
% Author: David S. White 
% Updated: 2021-10-08
% License: MIT
% ----------------------------------------------

% remove zero pairs
idxX = find(x==0);
idxY = find(y==0); 
idxZ = union(idxX, idxY); 

x(idxZ) = []; 
y(idxZ) = []; 

% nx = calcoptimalhistbins(x', 20, round(length(unique(x))));
% [~,xbins] = histcounts(x, nx); 
% bwX = log10(xbins(2)-xbins(1)); 
% 
% ny = calcoptimalhistbins(x', 20, round(length(unique(x))));
% [~,xbins] = histcounts(x, nx); 
% bwX = log10(xbins(2)-xbins(1)); 

% bw = 0.5;
[r,p] = corrcoef(x,y);
r = r(1,2); % simplify
p = p(1,2);
% minB = -1;
% maxB = 4; 
[counts, bins] = hist3(log10([x,y]), {minB:binWidth:maxB, minB:binWidth:maxB});
xx = min(bins{1}):binWidth:max(bins{1});
yy = min(bins{2}):binWidth:max(bins{2});
[X1, Y1] = meshgrid(bins{1}, bins{2});
[X2, Y2] = meshgrid(xx,yy); 
counts = interp2(X1, Y1, counts', X2, Y2); 

h = figure; 
contour(xx,yy,counts); 
colormap(jet); 
colorbar

if showPoints
    hold on
    scatter(log10(x), log10(y), 3,'filled', 'o', 'MarkerEdgeColor','k', 'MarkerFaceColor', 'k');
end

hold on 
plot(minB:1:maxB, minB:1:maxB, '--k')

% titleText = ['N = ', num2str(length(x)), ', {\it r} = ', num2str(round(r,2))];
% title(titleText)
%text(maxB-1.25, maxB-0.25, [['N = ', num2str(length(x))], newline,...
%    ['{\itr} = ', num2str(round(r,2))]], 'Fontsize',5)
% if p < 0.001
    text(minB+1, maxB+0.3, [['N = ', num2str(length(x))],...
        [', {\itr} = ', num2str(round(r,3))], ' ({\itp} = ', sprintf('%.0d',p), ')'], 'Fontsize',5)
% else
%     text(minB+0.2, maxB+0.2, [['N = ', num2str(length(x))],...
%         [', {\itr} = ', num2str(round(r,3))], ' ({\itp} = ', num2str(round(p, 3)), ')'], 'Fontsize',5)
% end
xticks(minB:1:maxB);
yticks(minB:1:maxB);
xlim([minB, maxB]);
ylim([minB, maxB]);

end