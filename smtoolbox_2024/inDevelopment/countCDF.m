function [x_unique,y_count,y_cdf] = countCDF(x)
% David S. White 
% 2023-08-01 
% quick function to get cummulative counts (c) and cummulative probability
% (p) from x [N,1]

x_unique = unique(x);
y_count = [cumsum(histcounts(x, x_unique)), numel(x)];
y_cdf = y_count./ numel(x);

end