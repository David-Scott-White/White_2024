function [dwellFits, yPDF] = fitDwells(dwells, nExp, bootStrap, fitLimits, tauGuess)
% David S. White
% dwhite7@wisc.edu
%
% Fit dwell time distributions with exponential(s). Merges and replaces
% previous fitExpDist.m and fitExpDist2.m from 01/21
%
% see plotDwells.m
%0

% Updates:
% --------
% 2021-02-04    DSW wrote code
%
% input:
% ------
% dwells = [Nx1] or [1xN] vector of dwell times
% nExp = max number of exponents to fit. If nExp == 3, fits 1 & 2.
%        default: nExp = 2; (Currently only fits 2 & 3 since 3 is unstable)
% nBoot = bootstrap fit. Will report 95% CI across nBoot in dwellFits
%           default: nBoot = 0 (does not perform bootstrapping)
%           note: diminishing returns nBoot > 1000
% fitLimits = [minDwellTime (frame rate), maxDwellTime (duration)]
%           default: fitLimits = [];
%     * see Kaur et al., 2019 Methods *
%
% output:
% ------
% dwellFits = struct with fits. Fields clearly named (IMO...)
% -------------------------------------------------------------------------

% check input, assign default
dwells = dwells(:);
nDwells = length(dwells);
if nDwells < 1
    dwellFits =[]; 
    yPDF = [];
    return 
end
if ~exist('nExp','var') || isempty(nExp); nExp = 2; end
if ~exist('bootStrap','var') || isempty(bootStrap); bootStrap = 0; nBoot=0; else; nBoot=1000; end
if ~exist('fitLimits', 'var') || isempty(fitLimits); fitLimits = []; end
if ~exist('tauGuess', 'var') || isempty(tauGuess); tauGuess = []; end

% output variable
dwellFits = struct;
dwellFits.nDwells = nDwells;
dwellFits.nBoot = nBoot;
dwellFits.limits = fitLimits;
dwellFits.tauGuess = tauGuess;

expPDF1 = @(t,tau1) 1/tau1*exp(-t/tau1);
expPDF2 = @(t,A1,tau1,tau2) A1.*expPDF1(t,tau1) + (1-A1).*expPDF1(t,tau2);

% set upper bound of fit
xmu = mean(dwells);
ub = max(dwells)*10;
lb = min(dwells)/10;

guessFactor = 2; 

% Fit parameters, bootstrap data?
phat = cell(nExp,1);
pci = cell(nExp,1);

% monoExp
[phat{1}, pci{1}] = mle(dwells,'pdf',expPDF1,'start', xmu,'lowerbound',lb,'upperbound',ub, ...
           'Options',statset('MaxIter',1e5, 'MaxFunEvals', 1e5));

% biExp
if nExp > 1
    if isempty(tauGuess)
        start = [0.5, mean(dwells)/guessFactor, mean(dwells)*guessFactor];
    else
        start = tauGuess;
    end
    [phat{2},pci{2}] = mle(dwells,'pdf',expPDF2,'start', start,...
        'lowerbound',[1e-3, lb, lb],'upperbound',[0.999, ub, ub], ...
           'Options',statset('MaxIter',1e5, 'MaxFunEvals', 1e5));
end

% constrain by fit limits? These equations are more sensitive to
% starting so run normal first, then constrain.
if exist('fitLimits', 'var') && ~isempty(fitLimits)
    tmin = fitLimits(1);
    tmax = fitLimits(2);
    expPDF1 = @(t,tau1) 1/(exp(-tmin/tau1)-exp(-tmax/tau1)) * (exp(-t/tau1)/tau1);
    expPDF2 = @(t, A1, tau1, tau2) (A1/tau1*exp(-t/tau1) + (1-A1)/tau2*exp(-t/tau2))/...
        ((A1*(exp(-tmin/tau1)-exp(-tmax/tau1))) + ((1-A1)*(exp(-tmin/tau2)-exp(-tmax/tau2))));
    
    % have to constrain the fit form emperical tests. Essentially, if we
    % assume that we do not have a conditional probabilty, fits should
    % always be an under estimate of true rate; therefore min values =
    % previous fit
    
    % mono
    [phat{1}, pci{1}] = mle(dwells,'pdf',expPDF1,'start', mean(dwells),'lowerbound',lb,'upperbound',ub, ...
           'Options',statset('MaxIter',1e5, 'MaxFunEvals', 1e5));
    
    % biExp
    if nExp > 1
        [phat{2}, pci{2}] = mle(dwells,'pdf',expPDF2,'start', [phat{2}(1), phat{2}(2), phat{2}(3)],...
           'lowerbound',[1e-3, lb, lb],'upperbound',[0.999, ub, ub],...
           'Options',statset('MaxIter',1e5, 'MaxFunEvals', 1e5));
    end
end

% NEED TO CORRECT THIS --------------------------------------------------
nBoot = 0;
if bootStrap
    nBoot = 1000;
    rng(42, 'twister'); % <- reproducibility
    idx = randi(length(dwells), length(dwells),nBoot,1); % may have error from memory setting...
    % monoExp
    phatTemp1 = zeros(nBoot,1);
    wb = waitbar(0, 'Bootstrapping Monoexponential...');
    for i = 1:nBoot
        phatTemp1(i,1) = mle(dwells(idx(:,i)),'pdf',expPDF1,'start', xmu,'lowerbound',lb,'upperbound',ub);
        waitbar(i/nBoot, wb);
    end
    dwellFits.bs{1} = phatTemp1;
    close(wb);
    
    % put in fix for if nboot == 1000
    phat{1} = mean(phatTemp1);
    if nBoot == 1000  % percentile bootstrap CI from boot strap
        [~,s] = sort(phatTemp1(:,1));
        pci{1} = [phatTemp1(s(25),:); phatTemp1(s(975),:)];
    else % bootstrap
        ci =  1.96 * std(phatTemp1)/sqrt(nDwells); % std to CI conversion
        pci{1} = [mean(phatTemp1)-ci; mean(phatTemp1)+ci];
    end
    
    % biExp
    if nExp > 1
        phatTemp2 = zeros(nBoot,3);
        % start = phat{2};
        wb = waitbar(0, 'Bootstrapping Biexponential...');
        for i= 1:nBoot
            x = dwells(idx(:,i));
            if isempty(tauGuess)
                start = [0.5, mean(x)/guessFactor, mean(x)*guessFactor];
            else
                start = tauGuess;
            end
            if start(2) < lb; start(2) = lb; end
            if start(3) > ub; start(3) = ub; end
            phatTemp2(i,1:3) = mle(dwells(idx(:,i)),'pdf',expPDF2,'start', [start(1), start(2), start(3)],...
                'lowerbound',[1e-3, lb, lb],'upperbound',[0.999, ub, ub],...
                'Options',statset('MaxIter',10000, 'MaxFunEvals', 10000));
            waitbar(i/nBoot, wb);
        end
        dwellFits.bs{2} = phatTemp2;
        close(wb);
        phat{2} = zeros(1,3);
        pci{2} = zeros(2,3);
        phatMu = mean(phatTemp2);
        phat{2} = phatMu;
        phatSd = std(phatTemp2);
        if nBoot == 1000 
            [~,s] = sort(phatTemp2(:,1));
            pci{2} = [phatTemp2(s(25),:); phatTemp2(s(975),:)];
            
        else
            ci = 1.96 * phatSd/(nDwells^0.5);
            pci{2}(1,:) = phatMu - ci;
            pci{2}(2,:) = phatMu + ci;
        end
    end
end

% Store values in output structure, compute LL of each fit
dwellFits.expPDF1 = expPDF1;
dwellFits.expPDF2 = expPDF2;
dwellFits.monoExpTau = phat{1};
if bootStrap
    dwellFits.monoExpSE = std(phatTemp1);
else
    dwellFits.monoExpSE = (phat{1} -pci{1}(1))/1.96; % estimate
end
dwellFits.monoExpCI = pci{1};
P1 = expPDF1(dwells, dwellFits.monoExpTau);
yPDF = P1;
P1(P1==0)=[];
dwellFits.monoExpPDF = P1;
dwellFits.monoExpLL = sum(log(P1));

if nExp > 1 % (sort with larger amplitude first)
    % if phat{2}(1) >= 0.5
    if phat{2}(2) < phat{2}(3)
        dwellFits.biExpAmp = [phat{2}(1),1-phat{2}(1)];
        dwellFits.biExpTau = [phat{2}(2),phat{2}(3)];
        if bootStrap
        dwellFits.biExpTauSE = [phatSd(2), phatSd(3)];
        dwellFits.biExpAmpSE = [phatSd(1), std(1-phatTemp2(:,1))];
        end
        dwellFits.biExpAmpCI = [pci{2}(:,1), flip(1-pci{2}(:,1))];
        dwellFits.biExpTauCI = [pci{2}(:,2), pci{2}(:,3)];  
        if ~bootStrap
            dwellFits.biExpAmpSE = (dwellFits.biExpAmp-dwellFits.biExpAmpCI(1,:))/1.96; 
            dwellFits.biExpTauSE = (dwellFits.biExpTau-dwellFits.biExpTauCI(1,:))/1.96;
        end 
    else
        dwellFits.biExpAmp = [1-phat{2}(1), phat{2}(1)];
        dwellFits.biExpTau = [phat{2}(3),phat{2}(2)];
        if bootStrap
            dwellFits.biExpTauSE = [phatSd(3), phatSd(2)];
            dwellFits.biExpAmpSE = [std(1-phatTemp2(:,1)),phatSd(1),];
        end
        dwellFits.biExpAmpCI = [flip(1-pci{2}(:,1)), pci{2}(:,1)];
        dwellFits.biExpTauCI = [pci{2}(:,3), pci{2}(:,2)];
        if ~bootStrap
            dwellFits.biExpAmpSE = (dwellFits.biExpAmp-dwellFits.biExpAmpCI(1,:))/1.96; 
            dwellFits.biExpTauSE = (dwellFits.biExpTau-dwellFits.biExpTauCI(1,:))/1.96;
        end
    end
    P2 = expPDF2(dwells,dwellFits.biExpAmp(1),dwellFits.biExpTau(1),dwellFits.biExpTau(2));
    yPDF=[yPDF,P2];
    P2(P2==0)=[];
    dwellFits.biExpPDF = P2;
    dwellFits.biExpLL = sum(log(P2));
    
    % Compare the fits.
    % 1. Loglikelihood ratio test (95% CI)
    
    % bootstreapped o
    LLR = 2*(dwellFits.biExpLL-dwellFits.monoExpLL);
    p = chi2cdf(LLR, 2, 'upper');
    %p = 1 - chi2cdf(LLR, 2);
    if p > 0.05
        best = 1;
    else
        best = 2;
    end
    dwellFits.LLR = [LLR, p, best];
    
    % 2. compute BIC
    BIC1 = -2*dwellFits.monoExpLL + log(nDwells);
    BIC2 = -2*dwellFits.biExpLL + 3*log(nDwells);
    if BIC2>=BIC1
        best = 1;
    else
        best = 2;
    end
    dwellFits.BIC = [BIC1, BIC2, best];
    
    % 3. Kullback-Leibler Divergence
    % data should be discrete bins, else use auto. Uses log instead of log2
    if length(dwells)>1
        dwellSorted = unique(sort(dwells));
        if dwellSorted(2) == dwellSorted(1)*2
            bins = dwellSorted(1):dwellSorted(1):dwellSorted(end);
            counts0 = hist(dwells, bins);
        else
            [counts0,bins] = histcounts(dwells,'BinMethod','auto');
            bins(end)=[];
        end
        areaAdjust = (bins(2)-bins(1))*sum(counts0);
        counts1 = expPDF1(bins, dwellFits.monoExpTau) * areaAdjust;
        counts2 = expPDF2(bins, dwellFits.biExpAmp(1), dwellFits.biExpTau(1), dwellFits.biExpTau(2)) * areaAdjust;
        
        % Probabilites (normalized values)
        p0 = counts0./(sum(counts0));
        p1 = counts1./(sum(counts1));
        p2 = counts2./(sum(counts2));
        
        % Fit w/ monoExp
        kl1 = p0 .* log(p0./p1);
        kl1(isnan(kl1))=0;
        kl1 = sum(kl1);
        
        % Fit w/ Bi Exp
        kl2 = p0 .* log(p0./p2);
        kl2(isnan(kl2))=0;
        kl2 = sum(kl2);
        
        if kl2>=kl1
            best = 1;
        else
            best = 2;
        end
        dwellFits.KLD = [kl1, kl2, best];
    else
        dwellFits.KLD = [];
    end
end

end
