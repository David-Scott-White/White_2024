function Kd = computeKd(deltaG, T)

% returns Kd in M

% Default temperature_k at Room temp if not proivded
if ~exist('T', 'var') || isempty(T)
    T = 293;
end
R = 1.9858775e-3; % gas constant in kcal/mol/K
Kd = exp(deltaG/(R*T)) ;


end