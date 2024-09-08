function y = geoSeries(x, p, n)
% Generate a geometric series
%
% Input : 
% x = first value 
% p = multiplicity factor 
% n = length of array 

% Output 
% y = vector of series
% 
% notes, probably an easier way to do this..
% 
% David S. White 
% 2022-01-06
% MIT
% -------------------------------------------------------------------------
if nargin < 3
    n = 50;
end
y = zeros(1,n); 
y(1) = x; 
for i = 2:n
    y(i) = y(i-1)*p;
end