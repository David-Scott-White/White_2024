function value_counts = countUnique(x, y)
% -------------------------------------------------------------------------
% Returns the unique values of x and associated counts
% -------------------------------------------------------------------------
%
% Adapted from https://www.mathworks.com/help/matlab/ref/double.unique.html
%
% Input:
%   x = array of values [N,1]
%   y = (Optional), return counts associated with array y [N,1]
%
% Output
%   value_counts = [values, counts];
%
% David S. White
% 2021-11-17
% MIT
% -------------------------------------------------------------------------

x = x(:);
[x_unique,~,ic] = unique(x);
x_counts = accumarray(ic,1);

if exist('y', 'var') && ~isempty(y)
    y = y(:);
    x_temp = zeros(length(y),1);
    p = 1; 
    for i = 1:length(y)
        if sum(y(i) == x_unique)
            x_temp(i) = x_counts(p);
            p = p +1;
        end
    end
    value_counts = [y, x_temp];
    
else
    value_counts = [x_unique, x_counts];
end
end