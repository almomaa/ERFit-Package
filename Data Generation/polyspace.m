function [A,P] = polyspace(X, order)

% Power Polynomials Expansion
%
%   This function is part of the Entropic Regression package (erfit)
%
%
% Inputs:
%    X: is an n-by-d matrix. n points in the d-dimensional sample space.
%    order: The expansion order.


% Outputs:
%    A: The expansion matrix (functions Library, dictionary,...)
%    P: The power matrix used to create A. where the ith row (P_i) is
%       d-dimensional vector, and the column A_i = Product_{j=1}^{d}
%       X{:,i}^P_{ij}.
%
%
% ---------------------------------------------------
%   Abd AlRahman R. AlMomani
%   aaalmoma(at)clarkson(dot)edu
%   April 2020
% ---------------------------------------------------

% For expansion of multivariate polynomials, the upper bound on the computing 
% time of the binary expansion is much greater than the upper bound on the 
% computing time of the iterative multiplication. See: LEE E. HEINDEL, 
% Computation of Powers of Multivariate Polynomials over the Integers, 
% JOURNAL OF COMPUTER AND SYSTEM SCIENCES 6, 1--8 (1971) 
% 
% So, it is more efficient to use iterative multiplication. However, it
% does not convienient to repeat cascading for loops for each expansion
% order.
%
% Here, we implemented the iterative multiplication up to the 6th expansion
% order. Higher expansion order can also be done for low dimension systems.

N = size(X,2);


if (order <= 6)
    [A, P] = lowOrder(X,order);
elseif ((N<=5) && (order<=10)) || (N<=3)
    P = nchoosek(kron(ones(1,N),0:order),N);
    P = unique(P,'rows');
    P(sum(P,2)>order,:) = [];
    P = flip(P,2);
    A = x2fx(X,P);
else
    msg = ['You request an expensive operation. ',...
           'For high dimension systems, the function return the',...
           ' expansion including only the terms up to the 6th order'];
    warning(msg)
    [A, P] = lowOrder(X,order);
end


end


function [A, P] = lowOrder(X,order)

[m,n] = size(X);
A = ones(m,1);
A = cat(2,A,X); % Add the constant and the linear terms
P = [zeros(1,n); eye(n)]; % Add power string of the constant 
                          % and the linear terms

if (order > 1), [A,P] = add2ndOrderTerms(A,P,X); end
if (order > 2), [A,P] = add3rdOrderTerms(A,P,X); end
if (order > 3), [A,P] = add4thOrderTerms(A,P,X); end
if (order > 4), [A,P] = add5thOrderTerms(A,P,X); end
if (order > 5), [A,P] = add6thOrderTerms(A,P,X); end


end

function [A, P] = add2ndOrderTerms(A,P,X)
N = size(X,2);
for i=1:N
    for j=i:N
        A = cat(2,A,prod([X(:,i) X(:,j)],2));
        ti = zeros(1,N); ti(i) = 1;
        tj = zeros(1,N); tj(j) = 1;
        t = ti + tj;
        P = cat(1,P,t);
        
    end
end
end

function [A,P] = add3rdOrderTerms(A,P,X)
N = size(X,2);
for i=1:N
    for j=i:N
        for k=j:N
            A = cat(2,A,prod([X(:,i) X(:,j) X(:,k)],2));
            ti = zeros(1,N); ti(i) = 1;
            tj = zeros(1,N); tj(j) = 1;
            tk = zeros(1,N); tk(k) = 1;
            t = ti + tj + tk;
            P = cat(1,P,t);
        end
    end
end
end

function [A,P] = add4thOrderTerms(A,P,X)
N = size(X,2);
for i=1:N
    for j=i:N
        for k=j:N
            for u=k:N
                A = cat(2,A,prod([X(:,i) X(:,j) X(:,k) X(:,u)],2));
                ti = zeros(1,N); ti(i) = 1;
                tj = zeros(1,N); tj(j) = 1;
                tk = zeros(1,N); tk(k) = 1;
                tu = zeros(1,N); tu(u) = 1;
                t = ti + tj + tk + tu;
                P = cat(1,P,t);
            end
        end
    end
end
end

function [A,P] = add5thOrderTerms(A,P,X)
N = size(X,2);
for i=1:N
    for j=i:N
        for k=j:N
            for u=k:N
                for l=u:N
                A = cat(2,A,prod([X(:,i) X(:,j) X(:,k) X(:,u) X(:,l)],2));
                ti = zeros(1,N); ti(i) = 1;
                tj = zeros(1,N); tj(j) = 1;
                tk = zeros(1,N); tk(k) = 1;
                tu = zeros(1,N); tu(u) = 1;
                tl = zeros(1,N); tl(l) = 1;
                t = ti + tj + tk + tu + tl;
                P = cat(1,P,t);
                end
            end
        end
    end
end
end

function [A,P] = add6thOrderTerms(A,P,X)
N = size(X,2);
for i=1:N
    for j=i:N
        for k=j:N
            for u=k:N
                for l=u:N
                    for n=l:N
                A = cat(2,A,prod([X(:,i) X(:,j) X(:,k) X(:,u) X(:,l) X(:,n)],2));
                ti = zeros(1,N); ti(i) = 1;
                tj = zeros(1,N); tj(j) = 1;
                tk = zeros(1,N); tk(k) = 1;
                tu = zeros(1,N); tu(u) = 1;
                tl = zeros(1,N); tl(l) = 1;
                tn = zeros(1,N); tn(n) = 1;
                t = ti + tj + tk + tu + tl + tn;
                P = cat(1,P,t);
                    end
                end
            end
        end
    end
end
end



