function [components, ideal, class] = runSKM(X, k)
% David S. White 
% 2022-11-2

% X = data
% k = number of states or starting guess (if length > 1)

tol = 1e-6;
maxIter = 100; 
done = 0; 
iter = 1;
llh = nan(maxIter,1);
while ~done
    if iter == 1
        [~, kmeans_class] = kmeansElkan(X, k(:));
        [viterbiPath, llh(iter)] = runViterbi(X, kmeans_class, 1);
        iter = iter + 1;
        
    else
        [viterbiPath, llh(iter)] = runViterbi(X, viterbiPath, 1);
        
        if ((llh(iter)) - (llh(iter-1))) <= tol || iter >= maxIter
            done = 1;
        else
            iter = iter + 1;
        end
    end
end
[components, ideal, class] = computeCenters(X, viterbiPath);
        
    