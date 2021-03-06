% Entropic Regression Backward Greedy Feature Elimination
%
%  This function is part of the Entropic Regression package (erfit)
%  
%
% Inputs:
%  A: Nxd basis matrix (Independent variables)
%  y: Nx1 dependent variable.
%  IX: the set of indices to investigate. (Elimination within these
%  indices)
%  options: Options structure. 
% 
%
% Outputs:
%  IX: The set of remaining indices after elimination.
%
%
% ---------------------------------------------------
%   Abd AlRahman R. AlMomani
%   aaalmoma(at)clarkson(dot)edu
%   May 2020
% ---------------------------------------------------
function IX = erbackward(A, y, IX, options)
val = -inf; ix = [];
keepin = options.keepin;

while ( (val <= options.tol) && (length(IX)>1) ) %Start the backward ER
    IX(ix) = [];           %Eliminate the selected weak candidate
    D = inf(1,length(IX)); %Initialize mutual information vector
    parfor (i=1:length(D),options.NumWorkers) %For the ith strong candidate
        if ~ismember(IX(i),keepin)
        rem = setdiff(IX,IX(i)); %#ok find all other candidates except i
        % Find the extra information added by the ith strong candidate 
        % given the all other strong candidates. (causation entropy). 
        f1 = A(:,IX)*pinv(A(:,IX))*y; %#ok
        f2 = A(:,rem)*pinv(A(:,rem))*y;
        D(i) = options.MIEstimator(f1, y, f2, options); %#ok 
        end
    end

    % Find the weakest candidate that have the minimum extra information
    [val,ix] = min(D);

    %If the minimum extra information (val) is less than the minimum
    %accepted influnce (tol)... then the function will be removed from the
    %strong candidates (see IX(ix) = []; above ), because that means 
    %the function is too weak and has no significant influnce compared to
    %the other candidates. Otherwise, that if val>tol, that mean the 
    %weakest candidate is strong, and have significant influnce, So,
    %terminate the backward step.
end
end