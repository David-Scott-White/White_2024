function [beta, plotData] = fitQuadraticBindingCurve(x, y, trace_M)
% Fit a binding curve from fluorescence (anisotropy, polarization etc..)
% data

% David S. White 
% 2022-08-18
% License: MIT

% x is ligand concentration in M
% y is output of binding (can be normalized, anisotropy, etc..) 
% trace_M = trace concnetration used for signal detection

% add in option for guesses *********

% references: 
% https://elifesciences.org/articles/57264
% https://stackoverflow.com/questions/66375656/non-linear-curve-fitting-fitnlm-and-lsqcurvefit-give-different-error-values
% https://www.urmc.rochester.edu/MediaLibraries/URMCMedia/labs/kielkopf-lab/documents/fluorescence_equations.pdf 

%% check data
y = y(:);
x = x(:);

% quadraticFxn = @(p, x) p(1)*(((x+trace_M+p(2)) - (sqrt((x+trace_M+p(2)).^2 - 4*c*x)))./(2*trace_M));
% float min and max values
quadraticFxn = @(p, x) p(1)+(p(2)-p(1))*(((x+trace_M+p(3)) - (sqrt((x+trace_M+p(3)).^2 - 4*trace_M*x)))./(2*trace_M));

% Perform the fit. 
% https://stackoverflow.com/questions/66375656/non-linear-curve-fitting-fitnlm-and-lsqcurvefit-give-different-error-values
opts = statset('TolFun',1e-10, 'robust', 'off');
nlm = fitnlm(x, y, quadraticFxn, [min(y), max(y), median(x)], 'Options', opts)

beta = nlm.Coefficients.Estimate;

% plot data
xplot = linspace(min(x), max(x), 1e4); 
yplot = quadraticFxn(beta, xplot);
plotData = [xplot(:), yplot(:)];

figure; 
plot(x, y, 'ko', xplot, yplot,'b-')
set(gca,'xscale', 'log');
xlabel('X [M]')
ylabel('Intensity')
title(['Kd (M): ' , num2str(beta(3))])

