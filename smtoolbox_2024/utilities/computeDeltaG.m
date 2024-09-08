function deltaG = computeDeltaG(forward_rate, reverse_rate, temperature_k)
% -------------------------------------------------------------------------
% Compute deltaG in (kcal/mol) for a transition 
% -------------------------------------------------------------------------
% 
% dG = -RT*ln(k1/ k2) 
%
% input: 
%   forward_rate. numeric. in 1/s
%   reverse_rate. numeric. in 1/s
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
R = 1.9858775e-3; % gas constant in kcal/mol/K
deltaG = -R * temperature_k * log(forward_rate./reverse_rate); 


