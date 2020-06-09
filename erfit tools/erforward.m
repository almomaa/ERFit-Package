% Entropic Regression Forward Greedy Feature Selection
%
%  This function is part of the Entropic Regression package (erfit)
%  
%
% Inputs:
%  A: Nxd basis matrix (Independent variables)
%  y: Nx1 dependent variable.
%  IXin: the set of indices to investigate. (Selection within these
%  indices)
%  options: Options structure. 
% 
%
% Outputs:
%  IX: The set of remaining indices after elimination.
%  success: Boolean flag indicate if the forward selection succeed.
%           If high degree of uncertinity detected, and too much parameters
%           selected, the success flag will be set to false. Then, the
%           backward elimination will be applied for all indices. 
%           
%
%
% ---------------------------------------------------
%   Abd AlRahman R. AlMomani
%   aaalmoma(at)clarkson(dot)edu
%   May 2020
% ---------------------------------------------------
function [IX, success] = erforward(A,y,IXin, options)
success = true;
ix = []; IX = unique([options.fkeepin(:); options.keepin(:)])';

N = size(A,2);
options.tol = tolEstimate(y,options);

Imax = options.MIEstimator(y,A*pinv(A)*y,options);

Done = false;
while ~Done %Start the forward ER
    IX = cat(2,IX,IXin(ix)); %Add the selected strong candidate. 
    Ilocal = options.MIEstimator(y,A(:,IX)*pinv(A(:,IX))*y);  
    
    D = -inf(1,N);      %Initialize mutual information vector
    parfor (i=1:N,options.NumWorkers)
        if ~ismember(IXin(i),IX) %Don't recheck what already selected
            % Find the extra information added by the ith candidate 
            % given the strong candidates selected so far
            f1 = A(:,[IX IXin(i)])*pinv(A(:,[IX IXin(i)]))*y; %#ok
            f2 = A(:,IX)*pinv(A(:,IX))*y;
            D(i) = options.MIEstimator(y,f1,f2,options); %#ok
        end
    end

    % Find the strongest candidate that have the maximum extra information
    [val,ix] = max(D);
    
    %If the maximum extra information (val) is more than the minimum   
    %accepted influnce (tol)... then the function will be added to the
    %strong candidates (see IX = [IX,IXin(ix)]; above ). Otherwise,
    %terminate the search.
    if ((Imax-Ilocal)<=options.tol) || (val<=options.tol)
        Done = true;
    elseif (length(IX)>max(8,N/2))
        success = false;
        Done = true;
    end



end

end