function [x_unique, x_counts] = countUniqueElements(x)

%https://www.mathworks.com/matlabcentral/answers/70061-how-to-count-unique-value-in-array

x = x(:);
x_unique = unique(x); 
x_counts = diff(find([true,diff(x')~=0,true]));
end