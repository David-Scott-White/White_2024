function [comps,ideal,class,metrics] = skm(x,varargin)
% David S. White
% dswhite2012@gmail.com

% updates:
% 2020-02-21    DSW wrote code
% 2021-02-23    DSW update with tolerance option

% overview:
% simple implementatin of the segmental k means algorithm for data
% idealization.
% This version features the ability to evaluate the likelihood of fitting
% with multiple states by use of an information theoretic criterion.

% example
% Q = [0,0.2, 0;0.2,0,0.2; 0,0.2,0];
% ss = simulateQ(Q,1:size(Q,1),10000,10);
% ssn = normrnd(ss,0.3);
% [comps,ideal,class, metrics] = skm(ssn,'min', 1,'max', 5,'iter',5,'fxn','BIC_GMM');

% check inputs
minK = 1;
maxK = 10;
maxIter = 100;
tol = 1e-6;
yStart = [];
oFxn = 'BIC_GMM';
for i = 1:2:length(varargin)-1
    switch varargin{i}
        case 'k'
            temp = varargin{i+1};
            if length(temp) == 1
                minK = temp(1);
                maxK = temp(1);
            else
                minK = temp(1);
                maxK = temp(2);
            end
        case 'start'
            yStart =  varargin{i+1};
            minK = length(unique(yStart));
            maxK = length(unique(yStart));
        case 'maxIter'
            maxIter = varargin{i+1};
        case 'fxn'
            oFxn = cell2mat(varargin(i+1));
        case 'tol'
            tol = varargin{i+1};
    end
end

% two options for how to optimize.
skmOption = 1;
if ~isempty(yStart)
    skmOption = 2;
end

% 1. Optimize across a given number of states

% perform clustering
nx = length(x);
totalK = maxK-minK+1;
kRange = minK:maxK;
viterbiPath = zeros(nx,totalK);
metrics = [];
switch skmOption
    case 1
        % optimze for minK to maxK, return the best across all values by
        % objective function
        for i = 1:totalK
            k = kRange(i);
            if k == 1
                viterbiPath = ones(nx,1);
            else
                llh = nan(maxIter,1);
                j = 1;
                done = 0;
                while ~done
                    if j == 1
                        [~,class] = kmeansElkan(x,k);
                        [vPath, llh(j)] = runViterbi(x,class, 1);
                        j = j+1;
                    else
                        comps = computeCenters(x, vPath);
                        [~,class] = kmeansElkan(x,comps(:,2));
                        [vPath, llh(j)] = runViterbi(x,class, 1);
                        if (llh(j) - llh(j-1)) <  tol || j == maxIter
                            done = 1;
                            viterbiPath(:,i) = vPath;
                        else
                            j = j +1;
                        end
                    end
                end
            end
        end
        % compute best fit
        [metrics, bestFit] = computeIC(x, viterbiPath, oFxn, 1);
        [comps,ideal,class] = computeCenters(x,viterbiPath(:,bestFit));
        
        % option 2, guess was provided.
    case 2
        % optimize using the provided start states (start can be array of
        % states or array of class assignments)
        if length(yStart) == length(x)
            j = 1;
            done = 0;
            llh = nan(maxIter,1);
            while ~done
                if j == 1
                    [vPath, llh(j)] = runViterbi(x, yStart, maxIter);
                    j = j +1;
                end
                comps = computeCenters(x, vPath);
                [~,class] = kmeansElkan(x,comps(:,2));
                [vPath, llh(j)] = runViterbi(x,class, 1);
                if j > 2
                    if (llh(j) - llh(j-1)) <  tol || j == maxIter
                        done = 1;
                    else
                        j = j +1;
                    end
                else
                    j = j +1;
                end
            end
        else
            % use provided values as start for kmeans
            j = 1;
            done = 0;
            llh = nan(maxIter,1);
            while ~done
                if j == 1
                    [~,class] = kmeansElkan(x,yStart(:));
                    [vPath, llh(j)] = runViterbi(x, class, maxIter);
                    j = j +1;
                end
                % need  way to check and store if fit improved..
                comps = computeCenters(x, vPath);
                [~,class] = kmeansElkan(x,comps(:,2));
                [vPath, llh(j)] = runViterbi(x, class, 1);
                if (llh(j) - llh(j-1)) <  tol || j == maxIter
                    done = 1;
                else
                    j = j +1;
                end
                
            end
        end
        figure; plot(llh)
        [comps,ideal,class] = computeCenters(x,vPath);
end
end
