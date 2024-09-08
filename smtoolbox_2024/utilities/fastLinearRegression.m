function [b, R2, RMSE, Yhat] = fastLinearRegression(x, y) 
% David S. White 
% 2021-10-01

% Applies the normal equation for linear regression
% input 
%  x = features [N,M]
%  y = observations [N,1]

X = [ones(length(x),1) x];
b = X\y;
Yhat = X*b; 
R2 = 1 - sum((y - Yhat).^2)/sum((y - mean(y)).^2);
RMSE = sqrt((sum((y - Yhat).^2))/length(y));
fprintf('Intercept:  %0.2d', b(1));
fprintf(' Slope:  %0.2d', b(2));
fprintf(' R2:  %0.2d', R2);
fprintf(' RMSE:  %0.2d', RMSE);
disp('')
end