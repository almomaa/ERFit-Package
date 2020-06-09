function fhandle = getODEHandle(P, S)
% Get ODE function handle
%
% Inputs:
%    P: Power Expansion Matrix
%    S: Sparse solution
 
% Outputs:
%    fhandle: ODE function handle
%
%
% ---------------------------------------------------
%   Abd AlRahman R. AlMomani
%   aaalmoma(at)clarkson(dot)edu
%   May 2020
% ---------------------------------------------------

F = cell(size(S,2),1);
for i=1:size(S,2)
    ix = ~~S(:,i);
    p = P(ix,:);
    B = S(ix,i);

    T = cell(1,size(p,1));
    for j=1:size(p,1)
        t = [];
        for k=1:size(p,2)
            if p(j,k)>1
                t = cat(2,t, '*x(', num2str(k), ')^',num2str(p(j,k)));
            elseif p(j,k)==1
                t = cat(2,t, '*x(', num2str(k),')');
            end
        end
        if isempty(t)
            t = num2str(B(j));
        elseif ((B(j)>0) && (j>1))
            t = cat(2,'+', num2str(B(j)), t);
        else
            t = cat(2,num2str(B(j)), t);
        end
        T{j} = t;
    end

    F{i} = ['@(x) ' cell2mat(T)];
end


F = cellfun(@str2func,F,'UniformOutput',false);
fhandle = @(t,x) cellfun(@feval,F,repmat({x},length(F),1));
