function [phat, pci, binomialCounts] = fitBinomial(observations, nTrials, zeroTruncated)
% David S. White
% dwhite7@wisc.edu
%
% Updates:
% 2020-04-23 DSW Wrote the code into function
% 2021-01-07 DSW. Comments and demo code.
%
% ---------
% Overview:
% ---------
% Maximum liklihood estimation of a binomial distribution for an array of
% integers. [Nx1] or [1xN]. Returns the probability of success (p) and
% expected distribution Common uses in single-molecule analysis
% include: 
%   1. Photobleaching step analysis (zeroTruncated = 1); 
%   2. State occupancy distributions (i.e., non-cooperative model). 
% 
% Requires Statistics and Machine Learning Toolbox 
%
% ---------
% Input:
% ---------
% observations = array [Nx1] or [1xN] of integers (0, 1, 2) 
% nTrials = integer. 
% zeroTruncated = boolean
% 
% ---------
% Output:
% ---------
% phat = estimated p of pdf (from mle) 
% pci = 95% confidence intervales of phat (from mle) 
% binomialCounts = [2,nTrials] of expected percents and counts

% ---------
% Demo Code: 
% ---------
% p = 0.3; 
% nTrials = 4; 
% zeroTruncated = 0; 
% observations = binornd(nTrials, p, 500,1);
% [phat, pci, binomialCounts] = fitBinomial(observations, nTrials, zeroTruncated);
% figure; bar(0:4, [accumarray(observations+1,1)'; binomialCounts(2,:)]); 
% xlabel('Counts'); ylabel('Observation'); legend('Observations', ['Fit: n=',num2str(phat(1)),'p=', num2str(phat(1))]);

%%  Check for zeroTruncated option
if nargin<3
    zeroTruncated = 0;
end

%% Evalute PDF
observations = observations(:);
binomialPDF = makeBinoPDF(nTrials, zeroTruncated); 
[phat,pci] = mle(observations, 'pdf', binomialPDF,'start',.5,'lowerbound',0,'upperbound',1);

%% make expexted counts
totalObservations = length(observations); 
binomialCounts = binomialPDF(zeroTruncated:nTrials,phat); % zeroTruncated = 0 or 1
binomialCounts = [binomialCounts; round(binomialCounts.*sum(totalObservations))];

end
%% Nested function to make pdf (not sure how else to change n value)
function binomialPDF = makeBinoPDF(nTrials, zeroTruncated)
if ~zeroTruncated
    if nTrials == 1; binomialPDF = @(x,p) binopdf(x,1,p);
    elseif nTrials == 2; binomialPDF = @(x,p) binopdf(x,2,p);
    elseif nTrials == 3; binomialPDF = @(x,p) binopdf(x,3,p);
    elseif nTrials == 4; binomialPDF = @(x,p) binopdf(x,4,p);
    elseif nTrials == 5; binomialPDF = @(x,p) binopdf(x,5,p);
    elseif nTrials == 6; binomialPDF = @(x,p) binopdf(x,6,p);
    elseif nTrials == 7; binomialPDF = @(x,p) binopdf(x,7,p);
    elseif nTrials == 8; binomialPDF = @(x,p) binopdf(x,8,p);
    elseif nTrials == 9; binomialPDF = @(x,p) binopdf(x,9,p);
    elseif nTrials == 10; binomialPDF = @(x,p) binopdf(x,10,p);
    end
else
    % make zero truncated binomiral
    if nTrials == 1; binomialPDF = @(x,p) binopdf(x,1,p) ./ (1-binopdf(0,1,p));
    elseif nTrials == 2; binomialPDF = @(x,p) binopdf(x,2,p) ./ (1-binopdf(0,2,p));
    elseif nTrials == 3; binomialPDF = @(x,p) binopdf(x,3,p) ./ (1-binopdf(0,3,p));
    elseif nTrials == 4; binomialPDF = @(x,p) binopdf(x,4,p) ./ (1-binopdf(0,4,p));
    elseif nTrials == 5; binomialPDF = @(x,p) binopdf(x,5,p) ./ (1-binopdf(0,5,p));
    elseif nTrials == 6; binomialPDF = @(x,p) binopdf(x,6,p) ./ (1-binopdf(0,6,p));
    elseif nTrials == 7; binomialPDF = @(x,p) binopdf(x,7,p) ./ (1-binopdf(0,7,p));
    elseif nTrials == 8; binomialPDF = @(x,p) binopdf(x,8,p) ./ (1-binopdf(0,8,p));
    elseif nTrials == 9; binomialPDF = @(x,p) binopdf(x,9,p) ./ (1-binopdf(0,9,p));
    elseif nTrials == 10; binomialPDF = @(x,p) binopdf(x,10,p) ./ (1-binopdf(0,10,p));
    end
end
end