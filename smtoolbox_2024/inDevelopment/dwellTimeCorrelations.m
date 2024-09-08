function [dwellPairs, r, n, p, xLabels, yLabels] = dwellTimeCorrelations(events, frame_rate_s)
% -------------------------------------------------------------------------
% Get Dwell Pairs and Correlate
% -------------------------------------------------------------------------
%
% input:
%   events = cell(N,1) from findEvents.m
%       currently only works for binary transitions of 1 <-> 2 (0 to 1?)
%
% output
%   dwellPairs = cell(N,1) where cell(i) = [M, 2] dwell pairs
%       % dwellPairs{1} = <bound> vs <unbound>
%
% Author: David S. White
% Updated: 2021-10-22
% License: MIT
% -------------------------------------------------------------------------

xLabels = {'log_1_0(\langleBound (s)\rangle)', 'log_1_0(Unbound_i (s))',...
    'log_1_0(Unbound_i (s))', 'log_1_0(Bound_i (s))'};
yLabels = {'log_1_0(\langleUnbound (s)\rangle)', 'log_1_0(Bound_i_+_1 (s))',...
    'log_1_0(Unbound_i_+_1 (s))', 'log_1_0(Bound_i_+_1 (s))'};

N = size(events, 1);
dwellPairs = cell(4,1);
r = zeros(1,4);
n = zeros(1,4);
p = zeros(1,4);

for i = 1:N
    e = events{i};
    nEvents = size(e,1);
    if nEvents > 1
        % Average Time per State ---
        unbound_mean = mean(e(e(:,4)==0, 3))*frame_rate_s;
        bound_mean = mean(e(e(:,4)==1, 3))*frame_rate_s;
        dwellPairs{1} = [dwellPairs{1}; [bound_mean, unbound_mean, i]];
    end
    % Unbound dwell vs bound dwell ---
    if nEvents > 1
        for j = 1:nEvents-1
            if e(j,4) == 0
                unbound_dwell = e(j,3)*frame_rate_s;
                bound_dwell = e(j+1,3)*frame_rate_s;
                dwellPairs{2} = [dwellPairs{2}; [unbound_dwell, bound_dwell]];
            end
        end
    end
    
    % unbound i vs unbound i+1
    if nEvents > 2
        for j = 1:nEvents-2
            if e(j,4) == 0
                if e(j+2,4) == 0
                    unbound_dwell_1 = e(j,3)*frame_rate_s;
                    unbound_dwell_2 = e(j+2,3)*frame_rate_s;
                    dwellPairs{3} = [dwellPairs{3}; [unbound_dwell_1, unbound_dwell_2]];
                end
            end
        end
    end
    
    % bound i vs bound i+1
    if nEvents > 2
        for j = 1:nEvents-2
            if e(j,4) == 1
                if e(j+2,4) == 1
                    bound_dwell_1 = e(j,3)*frame_rate_s;
                    bound_dwell_2 = e(j+2,3)*frame_rate_s;
                    dwellPairs{4} = [dwellPairs{4}; [bound_dwell_1, bound_dwell_2]];
                end
            end
        end
    end
end

% Get n and r for each dwell pair
for i = 1:4
    n(i) = size(dwellPairs{i},1); 
    [corr, prob] = corrcoef(dwellPairs{i});
    r(i) = corr(1,2); 
    p(i) = prob(1,2);
end

end
