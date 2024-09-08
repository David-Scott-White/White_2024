function [p0, rates] = evaluateQ(Q, stateLabels)
% -------------------------------------------------------------------------
% Evalaute the true Tau and Amplitudes of rates Q
% -------------------------------------------------------------------------
%
%
% 
% -------------------------------------------------------------------------

% the following is taken from simulateQ
% number of states;
K = size(Q,1);
states = unique(stateLabels);

if length(states) ~= K
    p0 = [];
    rates = []; 
    disp('Error in evaluateQ: Unique States and Number and Size of Q do not match')
    return
end

% unitary
q = Q; 
q = q-diag(sum(q,2));
S = [q ones(K,1)];

% equilibrium state probabilities
p0 = ones(1,K)/(S*(S'));

% invidual dwell times per state
ratesAll = zeros(K,1); 
for i = 1:K
    ratesAll(i) = 1/sum(Q(i,:))
end

end