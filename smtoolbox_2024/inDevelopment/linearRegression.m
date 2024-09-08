function [mdl, RMSE, R2] = linearRegression(x, y, varargin)
%% David S. White 
% 2023-09-08
% 
% Linear regression of x on y. Currently only 2D data
% requires matlab curve fitting toolbox 

% To plot results, see plotLinearRegression

% defaults;
slope = [];
intercept = [];
for i = 1:2:length(varargin)-1
    switch varargin{i}
        case 'slope'
            slope = varargin{i+1};
        case 'intercept'
            intercept = varargin{i+1};
    end
end

% ensure data is in the correct format
x = x(:);
y = y(:);

% normal equation to estimate parameters
X = [ones(length(x),1) x];
p0 = X\y;

% possible equations (of y = mx+b)

% fit both slope and intercept
if isempty(slope) && isempty(intercept)
    ft = fittype('m*x+b',...
        'dependent',{'y'},'independent',{'x'},...
        'coefficients',{'m','b'});
    mdl = fit(x, y, ft, 'StartPoint', [p0(2), p0(1)]);

% fit slope with provided intercept
elseif ~isempty(slope) && isempty(intercept)
    ft = fittype('m*x+b',...
        'dependent',{'y'},'independent',{'x'},...
        'coefficients', {'b'}, 'problem', 'm');
    mdl = fit(x, y, ft, 'StartPoint', p0(1), 'problem', slope);

% fit intercept with provided slope
elseif isempty(slope) && ~isempty(intercept)
     ft = fittype('m*x+b',...
        'dependent',{'y'},'independent',{'x'},...
        'coefficients',{'m'}, 'problem', 'b');
    mdl = fit(x, y, ft, 'StartPoint', p0(2), 'problem', intercept);

end

Yhat = feval(mdl, x); 
R2 = 1 - sum((y - Yhat).^2)/sum((y - mean(y)).^2);
RMSE = sqrt((sum((y - Yhat).^2))/length(y));

% print output
disp('Linear Regression:')
disp(['>> Slope: ', num2str(mdl.m)])
disp(['>> Intercept: ', num2str(mdl.b)])
disp(['>> R2: ', num2str(R2)])
disp(['>> RMSE: ', num2str(RMSE)])


end