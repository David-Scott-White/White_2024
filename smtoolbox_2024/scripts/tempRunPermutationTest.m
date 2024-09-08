% [p, observeddifference, effectsize] = permutationTest(sample1, sample2, permutations [, varargin])
N = 150; 
w = 0.9;
x = [exprnd(10, round(N*w),1); exprnd(200, round(N*(1-w)),1)];
y = [exprnd(10, round(N*w),1); exprnd(100, round(N*(1-w)),1)];
close all
permutationTest(x, y, 1e5, 'plotresult', 1, 'showprogress', 100);
