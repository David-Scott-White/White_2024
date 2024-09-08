function y = roundToSpecificNumber(x, n) 
% David S. White 
% 2022-08-05
% License: MIT 
%
%
% x = list of values 
% n = value to round to 
% 
% y = x rounded to nearest n 
%
% Round to nearest specific number 
% example: round to nearest "5"
% x = exprnd(100,100,1); 
% y = round( x / 5 ) * 5;
% [x,y]

y = round( x / n ) * n;
