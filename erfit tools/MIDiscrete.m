function I = MIDiscrete(x,y,varargin)

% Conditional Mutual Information estimator of a discrete random variable
%
% Inputs:
%    x: is an n-by-d1 symbolic matrix (that can be integers or characters)
%       n points in the d-dimensional sample space
%    y: is an n-by-d2 symbolic matrix (that can be integers or characters)
%       n points in the d-dimensional sample space
%    z: is an n-by-d3 symbolic matrix (that can be integers or characters)
%       n points in the d-dimensional sample space. (Optional Input) 
%
%    The function can be called with two variables (MIDiscrete(x,y)) to find the
%    mutual information. Providing the third input variable (MIDiscrete(x,y,z))
%    will find the conditional mutual information (the information between
%    x and y, given z).
 
% Outputs:
%   I = mutual information (if called with two variables).
%   I = conditional MI (if called with three variables).
%
% ---------------------------------------------------
%   Abd AlRahman R. AlMomani
%   aaalmoma(at)clarkson(dot)edu
%   April 2020
% ---------------------------------------------------

err = 'All imputs should have the same number of rows.';
assert((size(x,1)==size(y,1)), err);

if exist('eroptset','file')
    options = eroptset;
else
    options.BinMethod = auto;
end

if (~isempty(varargin)) && (isstruct(varargin{end}))
    options = varargin{end};
    varargin(end) = [];
end
Ni = length(varargin);
assert(ismember(Ni,[0 1]),'Error01: Number of inmputs')

N = size(x,1);


x = discretizeData(x,options);
y = discretizeData(y,options);

switch Ni
    case 0
        I = miDiscrete(x,y,options);
        
        if options.EmbeddedShuffleTest
            In = zeros(options.numPerm,1);
            for i=1:options.numPerm
                In(i) = miDiscrete(x,y(randperm(N)),options);
            end
            In = sort(In);
            tol = In(floor(options.alpha*options.numPerm));
            if I<tol, I = 0; end
        end
        
    case 1 
        z = varargin{1};
        z = discretizeData(z,options);
        assert((size(x,1)==size(z,1)), err);
        I = cmiDiscrete(x,y,z,options);
        
        if options.EmbeddedShuffleTest
            In = zeros(options.numPerm,1);
            for i=1:options.numPerm
                In(i) = cmiDiscrete(x,y(randperm(N)),z,options);
            end
            In = sort(In);
            tol = In(floor(options.alpha*options.numPerm));
            if I<tol, I = 0; end
        end
end

end

function I = miDiscrete(x,y,options)
    I = dEntropy(x) + dEntropy(y) - dEntropy([x,y]);
end

function I = cmiDiscrete(x,y,z,options)
    I = dEntropy([x,z]) + dEntropy([y,z]) - dEntropy(z) - ...
        dEntropy([x,y,z]);
end

