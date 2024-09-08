function [rate, rate_se] = tauToRate(tau, tau_se)
rate = 1./tau;
rate_se = rate - 1./(tau+tau_se);


