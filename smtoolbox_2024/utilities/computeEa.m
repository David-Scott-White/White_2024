function Ea = computeEa(k, A, temperature_k)
% -------------------------------------------------------------------------
% Compute deltaG in (kcal/mol) for a transition 
% -------------------------------------------------------------------------
% 
% k = Aexp(-Ea/RT)
%
% input: 
%   k. numeric. rate in 1/s
%   A. pre-exponetial term
%   temperature_k. numeric in K. default = 293 (room temp)
%
% output: 
%   deltaG. numeric. in kcal/mol
% 
% Author: David S. White 
% Updated: 2021-10-20
% License: MIT
% -------------------------------------------------------------------------

% Default temperature_k at Room temp if not proivded
if ~exist('temperature_k', 'var') || isempty(temperature_k)
    temperature_k = 293;
end
if ~exist('A', 'var') || isempty(A)
    temperature_k = 1;
end
R = 1.9858775e-3; % gas constant in kcal/mol/K
Ea = -log(k/A) * R * temperature_k;