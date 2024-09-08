function dG_ts = computeEyring(k, kappa, T_k)
% -------------------------------------------------------------------------
% Compute deltaG of transition state(dGt_ts) in (kcal/mol) for a transition 
% -------------------------------------------------------------------------
% 
% k = (kappa*kB*T)/h * exp(-dG_ts/RT)
%
% input: 
%   k. numeric. rate in 1/s
%   kappa. transmission coefficent. assumed unity (kappa = 1)
%   T_k. numeric in K. default = 293 (room temp)
%
% output: 
%   deltaG. numeric. in kcal/mol
% 
% Author: David S. White 
% Updated: 2021-10-20
% License: MIT
% -------------------------------------------------------------------------
if ~exist('T_k', 'var') || isempty(T_k)
    T_k = 293;
end
if ~exist('kappa', 'var') || isempty(kappa)
    kappa = 1;
end

R = 8.3145112;      % gas constant: 8.3145112 J K-1 mol-1 = 1.9872159 cal K-1 mol-1
kB = 1.380e-23;     % J/K(kB) = 3.297 × 10−24 cal/K = 1.380 × 10−23 J/K 
h = 6.6260755e-34;  % planck's constant (J s)

dG_ts_J = R*T_k*log(kappa*kB*T_k/h) - R*T_k * log(k);
dG_ts =  0.000239 * dG_ts_J; % kcal/mol output
