function [h, plot_data] = plotDwellsCDF(dwells, dwellFits, varargin)
% David S. White
% 2023-09-18

% takes results of fitDwells and plots
% see also plotDwells for histogram style plotting

% IN DEVELOPMENT... no options yets
markerSize = 2;

for i = 1:2:length(varargin)-1
    switch varargin{i}
        case 'markerSize'
            markerSize = varargin{i+1};
    end
end

dwells = dwells(:);
[x,~,y] = countCDF(dwells);
if isempty(dwellFits.limits)
    xx = linspace(x(1),x(end), 1e3);
    tau0 = dwellFits.monoExpTau;
    cdf1 = 1-exp(-xx/tau0);
    if isfield(dwellFits, 'biExpTau')
        amp1 = dwellFits.biExpAmp(1);
        tau1 = dwellFits.biExpTau(1);
        tau2 = dwellFits.biExpTau(2);
        cdf2 = 1-( (amp1*exp(-xx/tau1)) + ((1-amp1)*exp(-xx/tau2)) );
    else
        cdf2 = [];
    end
else
    xmin = dwellFits.limits(1);
    xmax = dwellFits.limits(2);
    % xx = xmin:xmin:xmax;
    xx = linspace(xmin, xmax, 1e3);
    tau0 = dwellFits.monoExpTau;
    k0 = 1/tau0;
    
    % cdf1 = 1-(-exp(-k0*(xx-xmax))/ (exp((xmax-xmin)*k0)-1));
    % cdf1 = cdf1-min(cdf1); % consider this the +C of the integral d
    
    
    cdf1 = (-(exp(-k0*xx)) / (exp(-k0*xmin)-exp(-k0*xmax)));
    %cdf1 = cdf1-min(cdf1);
    cdf1 = cdf1-cdf1(1);
    
    if isfield(dwellFits, 'biExpTau')
        amp1 = dwellFits.biExpAmp(1);
        amp2 = 1-amp1;
        tau1 = dwellFits.biExpTau(1);
        tau2 = dwellFits.biExpTau(2);
        k1 = 1/tau1;
        k2 = 1/tau2;
        %cdf2 = 1-(amp1*((exp(-k1*(xx-xmax))/ (exp((xmax-xmin)*k1)-1))) - (1-amp1)*((exp(-k2*(xx-xmax))/ (exp((xmax-xmin)*k2)-1))));
        %cdf2 = cdf2-min(cdf2);
        
        %cdf2 = ((amp1-1)*exp(-xx*k2) - (amp1)*exp(-xx*k1)) / (amp1*(exp(-k1*xmin)-exp(-k1*xmax)) + (1-amp1)*(exp(-k2*xmin)-exp(-k1*xmax)));
        %cdf2 = cdf2 - min(cdf2); % consider this the +C... the CDF is indefinate and does not behave well...
        
        cdf2 = (-amp1*exp(-xx*k1) - amp2*exp(-xx*k2)) / ( amp2*(exp(-k2*xmin)-exp(-k2*xmax))  + amp1*(exp(-k1*xmin)-exp(-k1*xmax)) );
        % cdf2 = cdf2-min(cdf2);
         cdf2 = cdf2-cdf2(1);
        
    else
        cdf2 = [];
    end
end

% Plot
h = figure;
hold on
scatter(x, y, markerSize, 'Marker', 'o', 'MarkerFaceColor',[0.5,0.5,0.5], 'MarkerEdgeColor', [0.5,0.5,0.5])
if ~isempty(cdf2)
    plot(xx, cdf1,'--', 'Color', [0 0.4470 0.7410])
    plot(xx, cdf2,'-', 'Color', 'r')
    
    nd = 1; %num decimals
    leg = legend(['N = ', num2str(length(dwells))],...
        ['\tau_0 = ', num2str(round(dwellFits.monoExpTau,nd)), ' ± ', num2str(round(dwellFits.monoExpSE,nd)), ' s'],...
        ['A_1 = ', num2str(round(dwellFits.biExpAmp(1),2)), ' ± ', num2str(round(dwellFits.biExpAmpSE(1),2)), ...
        newline, '\tau_1 = ', num2str(round(dwellFits.biExpTau(1),nd)), ' ± ', num2str(round(dwellFits.biExpTauSE(1),nd)), ' s', ...
        newline, '\tau_2 = ', num2str(round(dwellFits.biExpTau(2),nd)), ' ± ', num2str(round(dwellFits.biExpTauSE(2),nd)), ' s'],...
        'Location', 'northwest');
    leg.ItemTokenSize = leg.ItemTokenSize/3;
    
else
    plot(xx, cdf1,'-', 'Color', [0 0.4470 0.7410])
    nd = 1;
    leg = legend(['N = ', num2str(length(dwells))],...
        ['\tau_1 = ', num2str(round(dwellFits.monoExpTau,nd)), ' ± ', num2str(round(dwellFits.monoExpTau,nd)), ' s'],...
        'Location', 'southeast');
    leg.ItemTokenSize = leg.ItemTokenSize/3;
end
ylabel('Cumulative Probability')
xlabel('Time (s)')

plot_data=struct; 
plot_data.x = x(:); 
plot_data.y = y(:);
plot_data.xx = xx(:); 
plot_data.cdf1 = cdf1(:); 
plot_data.cdf2 = cdf2(:); 



