function [comps, ideal, class, metrics, all_ideal] = runPELT(X, varargin)
% David S. White
% 2022-11-02

% Run findchangepts from matlab with PELT algorithm. 
% Return in "DISC" format for ease. 
% Will also cluster data

% x = data; 

viterbi = 1; 
maxChangePoints = length(X)/2;
objectiveFunction = 'BIC_GMM';

for i = 1:2:length(varargin)-1
    switch varargin{i}
        case {'maxChangePoints', 'N'}
            maxChangePoints = varargin{i+1};
            
        case {'objectiveFunction', 'objFxn'}
            objectiveFunction = varargin{i+1};
        case 'viterbi'
            viterbi = varargin{i+1};
    end
end

ipt = findchangepts(X, 'MaxNumChanges',maxChangePoints,'Statistic', 'mean');
ipt =[1; ipt; length(X)];
ss = X*0;
x = length(ipt)-1;

for k = 1:length(ipt)-1
    s1 = ipt(k);
    s2 = ipt(k+1);
    ss(s1:s2) = x;
    x = x-1;
end

[~, ideal] = computeCenters(X, ss);
if ~isempty(objectiveFunction)
    all_ideal = aggCluster(X, ideal);
    [metrics, best_fit] = computeIC(X, all_ideal, objectiveFunction, 1);
    data_fit = all_ideal(:, best_fit);
else
    metrics = [];
    all_ideal = ideal; 
    data_fit = ideal;
end

[comps, ideal, class] = computeCenters(X, data_fit);

if viterbi
    data_fit  = runViterbi(X, class, 1);
    [comps, ideal, class] = computeCenters(X, data_fit);
end
