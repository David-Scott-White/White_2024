function x_scaled = scaleTimeSeries(x, x_fit)
% Scale raw time series of a.u. to [0,1,2,etc...]
% David S. White 
% 2020-12-28 

% input 
% x = raw time series [N,1] in a.u.
% x_fit = idealized values of x (i.e. disc_fit.ideal)

% note: 
% baseline state (min(x_fit)) is assumed to have a value of zero
[~,y_fit] = computeCenters(x, x_fit);
states = unique(y_fit);
nstates = length(states);
if nstates > 1
    x_scaled = (x-states(1))./(states(end)-states(1))*(nstates-1);
else
    x_scaled = (x./y_fit) -1 + x_fit; 
end

end