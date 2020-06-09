function S = lsfit(A,Y,ix)
% Ordinary Least Squares Solution
n = size(Y,2);
m = size(A,2);

S = zeros(m,n);
for i=1:n
    S(ix(:,i),i) = pinv(A(:,ix(:,i)))*Y(:,i);
end
end

