function E = dEntropy(X)

% Entropy estimator of a discrete random variable
%
% Inputs:
%    X: is an n-by-d symbolic matrix (that can be integers or characters)
%     n points in the d-dimensional sample space
 
% Outputs:
%    E: Entropy of X
%
%
% ---------------------------------------------------
%   Abd AlRahman R. AlMomani
%   aaalmoma(at)clarkson(dot)edu
%   May 2020
% ---------------------------------------------------

[~,~,ic] = unique(X,'rows');
P = accumarray(ic,1)./length(ic);
E = -dot(P,log2(P));
end

