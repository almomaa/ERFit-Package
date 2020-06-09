function [tol,In] = tolEstimate(y,options)

% Tolerence Estimation for mutual independence test.
%
%   This function is part of the Entropic Regression package (erfit)
%   
%
% Inputs:
%    y: is an n-by-d matrix. n points in the d-dimensional sample space.
%    options: Options structure that have the fields:
%      'numPerm': Number of permutation to perform in the permutation test.
%                 Default = 250.
%      'alpha': Confedence parameter. Default = 0.98.
%      'MIEstimator': Function handel for the mutual information estimator.
%
 
% Outputs:
%    tol: The tolerence value with confidence alpha.
%    In: numPerm x 1 vector that have all the values of the mutual
%    information for all permutations.
%
%
% ---------------------------------------------------
%   Abd AlRahman R. AlMomani
%   aaalmoma(at)clarkson(dot)edu
%   April 2020
% ---------------------------------------------------


N = size(y,1);
In = zeros(options.numPerm,1);
for i=1:options.numPerm
    In(i) = options.MIEstimator(y,y(randperm(N),:),options);
end
tol = quantile(In,options.alpha);
end

