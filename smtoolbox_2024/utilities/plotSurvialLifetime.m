function y = plotSurvialLifetime(x)
% David S. White 
% 2022-09-12
% 
%  x = data.rois structure anlyzed by smTraceViewer
%       Current assumotpions
%       > only looks at selected data
%       > expects selected and class == 1 being a molecule that never
%       photobleached / dissociated 
%
% 
% Returns y
%       [N,2] where N is number of frames
%       [N,1] = raw counts (ie molecules on)
%       [N,2] = normalized counts to N(1,1);

x = x(vertcat(x.status)==1);
n = length(x(1).fit.class);
y = zeros(n,2); % raw counts, normalized
for i = 1:length(x)
    if size(x(i).events,1) == 1
        % "unbound" state is actually bound
        y(:,1) = y(:,1) + (x(i).fit.class);
    else
        % subtract unbound state
         y(:,1) = y(:,1) + (x(i).fit.class)-1;
    end
end
y(:,2) = y(:,1)/ y(1,1);
% figure; plot(y(:,1))