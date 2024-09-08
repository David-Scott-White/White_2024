function [params, yhat, plotData] = fitBindingCurve(x, y, varargin)
% -------------------------------------------------------------------------
% Fit a binding curve 
% -------------------------------------------------------------------------
% Input
%   x = ligand concentration in M
%   y = fraction bound (in 0 to 1) 
%
% output
%   parms = struct
%       params.kd
%       params.bmax
%       params.hill
%       params.rmse
%   yhat = predicted fit values for provided xvalues
%
% David S. White 
% 2022-03-25 
% License: MIT 
% -------------------------------------------------------------------------

plotData = [];
yhat = []; 
params = []; 

% option argument check
for i = 1:2:length(varargin)-1
    switch varargin{i}
        case 'plot'
            plotData = varargin{i+1};
        case 'kd_M'
            kd_M = varargin{i+1};
        case 'bmax'
            bmax = varargin{i+1};
        case 'hill'
            hill = varargin{i+1};
        case 'constant'
            constant_M = varargin{i+1};
    end
end

% vary by input? 
% lets make a string for input type 
fxn = '111';

if ~exist('plotData', 'var'); plotData = 0; end
if ~exist('bmax', 'var'); bmax = []; fxn(1) = '0'; end
if ~exist('hill', 'var'); hill = [];fxn(2) = '0'; end
if ~exist('kd_M', 'var'); kd_M = []; fxn(3) = '0'; end
if exist('constant_M', 'var'); fxn = 'quadratic'; end

% if ~exist('constant_M', 'var'); constant_M = []; end

x = x(:); 
y = y(:); 
n = length(x); 

% log10 order of magnitude for x_range (e..g, 1e-6 returns -6);
nmin = floor(log10(min(x)));
nmax = ceil(log10(max(x)));

options = optimset('MaxIter', 1e5, 'MaxFunEvals', 1e5); 

% determine what to vary
switch fxn
    case '000'
        % constant: n/a
        % optimize: bmax, hill, kd
        bmax0 = max(y); 
        h0 = 1; 
        kd0 = median(x); 
        hillfxn = @(p, x) (p(1).*x.^p(2))./(x.^p(2)+p(3));
        costFunc = @(p) (sum((hillfxn(p, x) -  y).^2) / n).^0.5;
        p = fminsearch(costFunc, [bmax0, h0, kd0], options);
        bmax  = p(1); 
        hill = p(2); 
        kd_M = p(3); 
        yhat = (bmax.*x.^hill)./ (x.^hill+kd_M);
        
    case '100'
        % constant: bmax
        % optimize: hill, kd
        h0 = 1; 
        kd0 = median(x); 
        hillfxn = @(p, x, b) (b.*x.^p(1))./(x.^p(1)+p(2));
        costFunc = @(p) (sum((hillfxn(p, x, bmax) -  y).^2) / n).^0.5;
        p = fminsearch(costFunc, [h0, kd0], options);
        hill = p(1); 
        kd_M = p(2); 
        yhat = (bmax.*x.^hill)./ (x.^hill+kd_M);
        
        
    case '010'
        % constant: hill
        % optimize: bmax, kd
        bmax0 = max(y); 
        kd0 = median(x); 
        hillfxn = @(p, x, h) (p(1).*x.^h)./(x.^h+p(2));
        costFunc = @(p) (sum((hillfxn(p, x, hill) -  y).^2) / n).^0.5;
        p = fminsearch(costFunc, [bmax0, kd0], options);
        bmax = p(1); 
        kd_M = p(2);
        yhat = (bmax.*x.^hill)./ (x.^hill+kd_M);
        
    case '001'
        % constant: kd
        % optimize: bmax, hill
        bmax0 = max(y);
        h0 = 1;
        hillfxn = @(p, x, k) (p(1).*x.^p(2))./(x.^p(2)+k);
        costFunc = @(p) (sum((hillfxn(p, x, kd_M) -  y).^2) / n).^0.5;
        p = fminsearch(costFunc, [bmax0, h0], options);
        bmax = p(1);
        hill = p(2);
        yhat = (bmax.*x.^hill)./ (x.^hill+kd_M);
        
    case '110'
        % constant: bmax, hill
        % optimize: kd
        kd0 = median(x); 
        hillfxn = @(p, x, b, h) (b.*x.^h)./(x.^h+p(1));
        costFunc = @(p) (sum((hillfxn(p, x, bmax, hill) -  y).^2) / n).^0.5;
        p = fminsearch(costFunc, [kd0], options);
        kd_M = p(1);
        yhat = (bmax.*x.^hill)./ (x.^hill+kd_M); 
        
    case '101'
        % constant: bmax, kd
        % optimize: hill
        h0 = 1;
        hillfxn = @(p, x, b, k) (b.*x.^p(1))./(x.^p(1)+k);
        costFunc = @(p) (sum((hillfxn(p, x, bmax, kd_M) -  y).^2) / n).^0.5;
        p = fminsearch(costFunc, [h0], options);
        hill = p(1);
        yhat = (bmax.*x.^hill)./ (x.^hill+kd_M);
        
    case '011'
        % constant: hill, kd
        % optimize: bmax
        bmax0 = max(y);
        hillfxn = @(p, x, h, k) (p(1).*x.^h)./(x.^h+k);
        costFunc = @(p) (sum((hillfxn(p, x, hill, kd_M) -  y).^2) / n).^0.5;
        p = fminsearch(costFunc, [bmax0], options);
        bmax = p(1); 
        yhat = (bmax.*x.^hill)./ (x.^hill+kd_M);
         
    case '111'
        % constant: hill, bmax, kd
        % optimize: n/a (prediction only)
        yhat = (bmax.*x.^hill)./ (x.^hill+kd_M);
        
    case 'quadratic-bmax1'
        % assume hill and bmax = 1; (can vary bmax!)
        % optimize kd
        kd0 = median(x);
        quadraticFxn = @(p, x, c) (((x+c+p) - (sqrt((x+c+p).^2 - 4*c*x)))./(2*x));
        costFunc = @(p) (sum((quadraticFxn(p, x, constant_M) -  y).^2) / n).^0.5;
        kd_M = fminsearch(costFunc, kd0, options);
        xplot = logspace(nmin, nmax,100);
        xplot = xplot(:); 
        yhat = ((x+constant_M+kd_M) - (sqrt((x+constant_M+kd_M).^2 - 4*constant_M*x)))./(2*constant_M);
        yplot = ((xplot+constant_M+kd_M) - (sqrt((xplot+constant_M+kd_M).^2 - 4*constant_M*xplot)))./(2*constant_M);
        plotData = [xplot, yplot];
        
        
    case 'quadratic'
        % constant: constant_M, hill=1
        % optimze bmax and kd
        hill = 1; 
        bmax0 = max(y);
        kd0 = median(x);
        
        quadraticFxn = @(p, x, c) p(1)*(((x+c+p(2)) - (sqrt((x+c+p(2)).^2 - 4*c*x)))./(2*c));
        
        costFunc = @(p) (sum((quadraticFxn(p, x, constant_M) -  y).^2) / n).^0.5;
        p = fminsearch(costFunc, [bmax0, kd0], options);
        bmax = p(1);
        kd_M = p(2);
       
        xplot = logspace(nmin, nmax,100);
        xplot = xplot(:); 
        yhat = bmax*(((x+constant_M+kd_M) - (sqrt((x+constant_M+kd_M).^2 - 4*constant_M*x)))./(2*constant_M));
        yplot = bmax*(((xplot+constant_M+kd_M) - (sqrt((xplot+constant_M+kd_M).^2 - 4*constant_M*xplot)))./(2*constant_M));
        plotData = [xplot, yplot];
        
    case 'quadratic2'
        % constant: constant_M, hill=1
        % optimze bmax and kd
        hill = 1;
        bmax0 = max(y);
        kd0 = median(x);
        
        quadraticFxn = @(p, x, c) p(1)*(((x.^p(3)+c+p(2)) - (sqrt((x.^p(3)+c+p(2)).^2 - 4*c*x.^p(3))))./(2*c));
        
        costFunc = @(p) (sum((quadraticFxn(p, x, constant_M) -  y).^2) / n).^0.5;
        p = fminsearch(costFunc, [bmax0, kd0, 1], options);
        bmax = p(1);
        kd_M = p(2);
        hill = p(3);
        
        xplot = logspace(nmin, nmax,100);
        xplot = xplot(:);
        yhat = bmax*(((x+constant_M+kd_M) - (sqrt((x+constant_M.*+kd_M).^2 - 4*constant_M*x)))./(2*constant_M));
        yplot = bmax*(((xplot+constant_M+kd_M) - (sqrt((xplot+constant_M+kd_M).^2 - 4*constant_M*xplot)))./(2*constant_M));
        plotData = [xplot, yplot];
        
end
rmse = sqrt(sum((y-yhat).^2)/n); 

params = struct; 
params.kd_M = kd_M; 
params.bmax = bmax; 
params.hill = hill;
params.rmse = rmse; 




