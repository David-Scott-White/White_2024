function [mdl,yhat, plot_data] = fitExpCurve(x, y, t0, plotFigure)
% David S. White 
% 2022-09-27

% Constrain exp1 from fit(x,y,'exp1') to account for last value not being
% zero582318/////////////////
% Adapted from https://www.mathworks.com/matlabcentral/answers/733368-using-the-fit-function-to-find-an-exponential-gives-me-a-straight-line

x = x(:);
y = y(:);
if ~exist('t0', 'var') || isempty(t0)
    t0 = max(diff(y))/mean(diff(x));
end

if ~exist('plotFigure', 'var') || isempty(plotFigure)
    plotFigure = 0;
end

% free c
ft = fittype('a*exp(-b*t) + c','indep','t');
fitOptions = fitoptions('Method', 'NonlinearLeastSquares','StartPoint',[y(1), t0, max(y)]);



mdl = fit(x, y, ft, fitOptions);
yhat = mdl.a*exp(-mdl.b*x) + mdl.c;
xrange = linspace(x(1), x(end),1000);
yhat_plot = mdl.a*exp(-mdl.b*xrange) + mdl.c;

plot_data = [xrange; yhat_plot]';
if plotFigure
    h =figure;
    hold on
    plot(xrange,yhat_plot,'-r')
    plot(x,y,'k.');
    xlabel('x');
    ylabel('y')
    str1 = sprintf('b = %0.3d', mdl.b);
    str2 = sprintf('c = %0.2f', mdl.c);
    title([str1, ' | ', str2]);
end