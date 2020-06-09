function S = iterativeRWLS(A,Y,ix)
% Iterative Re-Weighted Least Squares Solution
n = size(Y,2);
m = size(A,2);

S = zeros(m,n);
for i=1:n
    S(ix(:,i),i) = robustfit(A(:,ix(:,i)),Y(:,i),'logistic',[],'off');
end
end



