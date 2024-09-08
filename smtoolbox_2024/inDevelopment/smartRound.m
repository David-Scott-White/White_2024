function xRound = smartRound(x)
% David S. White 
% 2021-12-29
% Rounds values based on order of magnitude

y = 0;
if x == 0
    disp('Error in smartround: x cannot be 0');
    return
end

x = abs(x);
n = floor(log10(x));
xRound = round(x, -1*n);

% xRound = ceil(x/10)*10;
end





