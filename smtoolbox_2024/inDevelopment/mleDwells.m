function dwellMLE = mleDwells(dwells, nExponentials, varargin)
% David S. White
% 2023-03-24
% Rewrite of former "fitDwells". Name changed for logic. Now up to 3 exp
%
%
% input:
%   dwells = dwells [N,1] or [1xN];
%   nExponentials = double. Accepted values = 1,2,3
%
% WORK IN PROGRESS

% -------------------------------------------------------------------------
dwells = dwells(:);
nDwells = length(dwells);

if ~exist('nExponentials', 'var') || isempty(nExponentials)
    nExponentials = 1;
elseif nExponentials > 3
    nExponentials = 3;
end

% Default Values for optional arguments -----------------------------------
bootstrap = 0; % boolean. 1000 bootstrap.
limits = []; % min Time, max time
paramstart = {};
maxIter = 1e5;
maxFunEvals = 1e5;
for i = 1:2:length(varargin)-1
    switch varargin{i}
        case 'bootstrap'; bootstrap = varargin{i+1};
        case 'dwellLimits'; limits = varargin{i+1};
        case 'start'; paramstart = varargin{i+1};
        case 'maxIter'; maxIter = varargin{i+1};
        case 'maxFunEvals'; maxFunEvals = maxFunEvals{i+1};
    end
end

% Initalize output stucture
dwellMLE = struct;
dwellMLE.bootstrap = bootstrap;
dwellMLE.limits = limits;
% parameters stored as cells for each fit
dwellMLE.pdf = cell(nExponentials,1);
dwellMLE.phat = cell(nExponentials,1); % parameter estimates
dwellMLE.pci = cell(nExponentials,1); % 95% confidence intervals
dwellMLE.se = cell(nExponentials,1);  % standard error of the parameters
dwellMLE.ll = cell(nExponentials,1);  % loglikelihood of the fit
if nExponentials > 1
    dwellMLE.LLR = cell(nExponentials-1,1);
    dwellMLE.BIC = cell(nExponentials-1,1);
else
    dwellMLE.LLR = [];
    dwellMLE.BIC = [];
end

% Define Functions
if isempty(limits)
    if nExponentials >= 1
        dwellMLE.pdf{1} = @(t,tau1) exppdf(t,tau1);
    end
    if nExponentials >= 2
        dwellMLE.pdf{2} = @(t,A1,tau1,tau2) A1.*exppdf(t,tau1) + (1-A1).*exppdf(t,tau2);
    end
    if nExponentials >= 3
        dwellMLE.pdf{3} = @(t,A1,A2,tau1,tau2,tau3) A1.*exppdf(t,tau1) + (A2).*exppdf(t,tau2) + (1-A1-A2).*exppdf(t,tau3);
    end
else
    % see Kaur et al., 2019 Methods
    tmin = limits(1);
    tmax = limits(2);
    if nExponentials >= 1
        dwellMLE.pdf{1} = @(t,tau1) 1/(exp(-tmin/tau1)-exp(-tmax/tau1)) * (exp(-t/tau1)/tau1);
    end
    if nExponentials >= 2
        dwellMLE.pdf{2} = @(t, A1, tau1, tau2) (A1/tau1*exp(-t/tau1) + (1-A1)/tau2*exp(-t/tau2))/...
            ((A1*(exp(-tmin/tau1)-exp(-tmax/tau1))) + ((1-A1)*(exp(-tmin/tau2)-exp(-tmax/tau2))));
    end
    if nExponentials >= 3
        dwellMLE.pdf{3} = @(t, A1, S2, tau1, tau2, tau) (A1/tau1*exp(-t/tau1) + (A2)/tau2*exp(-t/tau2) + (1-A1-A2)/tau3*exp(-t/tau3))/...
            ((A1*(exp(-tmin/tau1)-exp(-tmax/tau1))) + (A2*(exp(-tmin/tau2)-exp(-tmax/tau2))) + ((1-A1-A2)*(exp(-tmin/tau3)-exp(-tmax/tau3))));
    end
end

% Fit options
lowerboundAmp = 1e-3;
upperboundAmp = 1-1e-3;
lowerboundTau = min(dwells)/10;
upperboundTau = max(dwells)*10;
muDwells = mean(dwells);
options = statset('MaxIter',maxIter, 'MaxFunEvals', maxFunEvals);

if bootstrap
    rng(42, 'twister');                    % reproducibility
    idx = randi(nDwells, nDwells, 1000,1); % may have error from memory setting...
    dwells = dwells(idx);
end

% Maximum likelihood estimates --------------------------------------------
for i = 1:nExponentials

    % single exponential
    if i == 1
        p0 = muDwells;
        lb = lowerboundTau;
        ub = upperboundTau;
        if ~bootstrap
            [dwellMLE.phat{i}, dwellMLE.pci{i}] = mleExp1(dwells, dwellMLE.pdf{i}, p0, lb, ub, options);
            dwellMLE.se{i} = (dwellMLE.phat{i}-dwellMLE.phat{i}(1))/1.96;
        else
            phat = zeros(1000,3);
            for j = 1:1000
                phat(j,:) = mleExp1(dwells(:,j), dwellMLE.pdf{i}, p0, lb, ub, options);
            end
            [dwellMLE.phat{i}, dwellMLE.se{i}] = normfit(phat);
            phat = sort(phat);
            dwellMLE.pci{i} = [phat(25,:);phat(975,:)];
        end
        dwellMLE.ll{i} = evalExp1(dwells, dwellMLE.pdf{i}, dwellMLE.phat{i});
        
        % double exponential
    elseif i == 2
        p0 = [0.6, muDwells/3, muDwells*3];
        lb = [lowerboundAmp, lowerboundTau, lowerboundTau];
        ub = [upperboundAmp, upperboundTau, upperboundTau];
        if ~bootstrap
            [dwellMLE.phat{i}, dwellMLE.pci{i}] = mleExp2(dwells, dwellMLE.pdf{i}, p0, lb, ub, options);
            dwellMLE.se{i} = (dwellMLE.phat{i}-dwellMLE.phat{i}(1))/1.96;
        else
            phat = zeros(1000,3);
            for j = 1:1000
                phat(j,:) = mleExp2(dwells(:,j), dwellMLE.pdf{i}, p0, lb, ub, options);
            end
            [dwellMLE.phat{i}, dwellMLE.se{i}] = normfit(phat);
            phat = sort(phat);
            dwellMLE.pci{i} = [phat(25,:);phat(975,:)];
        end
        dwellMLE.ll{i} = evalExp2(dwells, dwellMLE.pdf{i}, dwellMLE.phat{i});
        
        % triple exponential
    elseif i == 3
        % temp 
        p0 = [0.2, 0.1, 300, 300, 100];
        lb = [lowerboundAmp, lowerboundAmp, lowerboundTau, lowerboundTau, lowerboundTau];
        ub = [upperboundAmp, upperboundAmp, upperboundTau, upperboundTau, upperboundTau];
        if ~bootstrap
            [dwellMLE.phat{i}, dwellMLE.pci{i}] = mleExp3(dwells, dwellMLE.pdf{i}, p0, lb, ub, options);
            dwellMLE.se{i} = (dwellMLE.phat{i}-dwellMLE.phat{i}(1))/1.96;
        else
            phat = zeros(1000,3);
            for j = 1:1000
                phat(j,:) = mleExp2(dwells(:,j), dwellMLE.pdf{i}, p0, lb, ub, options);
            end
            [dwellMLE.phat{i}, dwellMLE.se{i}] = normfit(phat);
            phat = sort(phat);
            dwellMLE.pci{i} = [phat(25,:);phat(975,:)];
        end
        dwellMLE.ll{i} = evalExp3(dwells, dwellMLE.pdf{i}, dwellMLE.phat{i});
    end
    
end

    
    
    % LLR & BIC
    % if nExponentials > 1
    %     LLR = zeros(nExponentials-1, 3);
    %     BIC = zeros(nExponentials-1, 3);
    %     nLogDwells = log(length(dwells));
    %     for i = 2:nExponentials
    %         k1 = length(phat{i-1});
    %         k2 = length(phat{i-1});
    %         % LLR
    %         llr = 2*(LL(i)-LL(i-1));
    %         p = 1 - chi2cdf(llr, k2-k1);
    %         if p > 0.05
    %             best = i -1;
    %         else
    %             best = i;
    %         end
    %         LLR(i-1,:) = [llr, p, best];
    %
    %         % BIC
    %         bic1 = -2*LL(i-1) + k1*nLogDwells;
    %         bic2 = -2*LL(i) + k2*nLogDwells;
    %         if bic2>=bic1
    %             best = i-1;
    %         else
    %             best = i;
    %         end
    %         BIC(i-1,:) = [bic1, bic2, best];
    %     end
    % else
    %     LLR = [];
    %     BIC = [];
    % end
    LLR = [];
    BIC = [];
    
    % compute BIC
    
    % add bootstrap options
    
    % sort output by taus (min to max)
    
    % evaluate the fits (LLH)
    
    % compare the fits (LLR, BIC)
    
end

% Local functions for maximum likelihood estimates ------------------------
    function [phat, pci] = mleExp1(x, pdfFunc, param0, lowerbound, upperbound, options)
        % phat = estimated paramters.
        % pci = confidence intervals
        
        [phat, pci] = mle(x,...
            'pdf', pdfFunc,...
            'start', param0,...
            'lowerbound', lowerbound,...
            'upperbound', upperbound, ...
            'Options', options);
        
    end

    function [phat, pci] = mleExp2(x, pdfFunc, param0, lowerbound, upperbound, options)
        % phat = estimated paramters.
        % pci = confidence intervals
        
        [phat, pci] = mle(x,...
            'pdf', pdfFunc,...
            'start', param0,...
            'lowerbound',[lowerbound(1),lowerbound(2),lowerbound(3)], ...
            'upperbound',[upperbound(1),upperbound(2),upperbound(3)], ...
            'Options', options);
        
        % sort phat and pci by tau (smallest first);
        if phat(2) > phat(3)
            phat = [1-phat(1), phat(3), phat(2)];
            pci = [pci(:,1), pci(:,3), pci(:,2)];
        end
    end
    function [phat, pci] = mleExp3(x, pdfFunc, param0, lowerbound, upperbound, options)
        % phat = estimated paramters.
        % pci = confidence intervals
        
        [phat, pci] = mle(x,...
            'pdf', pdfFunc,...
            'start', param0,...
            'lowerbound',[lowerbound(1),lowerbound(2),lowerbound(3),lowerbound(4),lowerbound(5)],...
            'upperbound',[upperbound(1),upperbound(2),upperbound(3),upperbound(4),upperbound(5)],...
            'Options', options);
        % sort phat and pci by tau (smallest first);
    end

% Compute loglikelihood of the fits -----------------------------------
    function ll = evalExp1(x, pdfFunc, params)
        y = pdfFunc(x, params(1));
        y(y==0)=[];
        ll = sum(log(y));
    end
    function ll = evalExp2(x, pdfFunc, params)
        y = pdfFunc(x, params(1), params(2), params(3));
        y(y==0)=[];
        ll = sum(log(y));
    end
    function ll = evalExp3(x, pdfFunc, params)
        y = pdfFunc(x, params(1), params(2), params(3), params(4), params(5));
        y(y==0)=[];
        ll = sum(log(y));
    end

