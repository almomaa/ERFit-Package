%% Conditional Mutual Information
% _cmiKnn_
% Finds the conditional mutual information $I(x;y|z)$ (reads, mutual
% information between $x$ and $y$ given $z$).

%% Citation
%  This code is part of the Entropic Regression Package.
%  This code is based of the theory described in several works.
%  If you use this code, consider giving the credit for all the work
%  that forms the theoritical base of this code.
%
%  To cite this code:
%
%  [ 1 ] Abd AlRahman R. AlMomani, Jie Sun, and Erik Bollt. How Entropic 
%        Regression Beats the Outliers Problem in Nonlinear System 
%        Identification. Chaos 30, 013107 (2020).
%
%  [ 2 ] Alexander Kraskov, Harald St¨ogbauer, and Peter Grassberger. 
%        Estimating mutual information. Physical Review E, 69:066-138,2004
%  
%  [ 3 ] Alexander Kraskov, Harald St¨ogbauer, and Peter Grassberger. 
%        Estimating mutual information. Physical Review E, 69:066-138,2004
%
%  [ 4 ] Alexander Kraskov, Harald St¨ogbauer, and Peter Grassberger. 
%        Estimating mutual information. Physical Review E, 69:066-138,2004
%

%% Inputs
%
% * x : Mxd column wise matrix (source).
% * y : Mxr column wise matrix (distination).
% * z : Mxk column wise matrix of condition set.
% * options: Structure that has two fields: 
%            - options.K: the integer value K for the Kth 
%              nearest neighbor to be used in estimation.
%            - options.distFun: cell array that have name of the 
%              distance function to be used, which should match the
%              built-in distance function defined by matlab (see
%              pdist2)
%

%% Outputs
% _I_ : scalar of the result conditioned mutual information




%% Function Body
%
%%
function [I] = MIKnn(x,y,varargin)

err = 'All imputs should have the same number of rows.';
assert((size(x,1)==size(y,1)), err);

if exist('eroptset','file')
    options = eroptset;
else
    options.K = 2;
    options.p = inf;
end

if (~isempty(varargin)) && (isstruct(varargin{end}))
    options = varargin{end};
    varargin(end) = [];
end
Ni = length(varargin);
assert(ismember(Ni,[0 1]),'Error01: Number of inmputs')

if (size(x,1)>2000) && options.EmbeddedShuffleTest
    warning(['Knn estimator with embedded shuffle test and large number '...
             'of measurements is expensive, and may take long time. ' ...
             'consider using discrete estimator instead'])
elseif size(x,1)>9000
    warning(['Knn Estimator with large number of measurements maybe ' ...
             'expensive, consider using discrete estimator instead'])
end


N = size(x,1);

switch Ni
    case 0
        I = miKnn(x,y,options);
        
        if options.EmbeddedShuffleTest
            In = zeros(options.numPerm,1);
            for i=1:options.numPerm
                In(i) = miKnn(x,y(randperm(N)),options);
            end
            In = sort(In);
            tol = In(floor(options.alpha*options.numPerm));
            if I<tol, I = 0; end
        end
        
    case 1 
        z = varargin{1};
        assert((size(x,1)==size(z,1)), err);
        I = cmiKnn(x,y,z,options);
        
        if options.EmbeddedShuffleTest
            In = zeros(options.numPerm,1);
            for i=1:options.numPerm
                In(i) = cmiKnn(x,y(randperm(N)),z,options);
            end
            In = sort(In);
            tol = In(floor(options.alpha*options.numPerm));
            if I<tol, I = 0; end
        end
end
end

function [I] = cmiKnn(x,y,z,options)

K = options.K;
distInfo = {'minkowski',options.p};


% If the condition set $z$ is empty, then use the Mutual inforation
% estimator.
if isempty(z)
    [I] = miKnn(x,y,options);
    return
end
 
% To construct the Joint Space between all variables we have:
JS = cat(2,x,y,z);
% Find the K^th smallest distance in the joint space JS = (x,y,z)
D = pdist2(JS,JS,distInfo{:},'Smallest',K+1)';
epsilon = D(:,end);
% Instead of the above two lines, the one may use the knnsearch function,
% but we found the above implementation is faster.

% Find number of points from $(x,z), (y,z)$, and $(z,z)$ that lies withing the
% K^{th} nearest neighbor distance from each point of themself.
Dxz = pdist2([x,z],[x,z],distInfo{:});
nxz = sum(bsxfun(@lt,Dxz,epsilon),2) - 1;

Dyz = pdist2([y,z],[y,z],distInfo{:});
nyz = sum(bsxfun(@lt,Dyz,epsilon),2) - 1;

Dz = pdist2(z,z,distInfo{:});
nz = sum(bsxfun(@lt,Dz,epsilon),2) - 1;

% VP Estimation formula.
I = psi(K) - mean(psi(nxz+1)+psi(nyz+1)-psi(nz+1)); 
end

function I = miKnn(x,y,options)
% Initialize and verify inputs.
K = options.K;
distInfo = {'minkowski',options.p};


% To construct the Joint Space between all variables we have:
JS = [x,y];  n = size(JS,1);


% Find the K^th smallest distance in the joint space JS
D = pdist2(JS,JS,distInfo{:},'Smallest',K+1)';
epsilon = D(:,end); %Set threshold value

% Find points on x with pairwise distance
% less than threshold value
Dx = pdist2(x,x,distInfo{:});
nx = sum(bsxfun(@lt,Dx,epsilon),2) - 1;

% Find points on y with pairwise distance
% less than threshold value
Dy = pdist2(y,y,distInfo{:});
ny = sum(bsxfun(@lt,Dy,epsilon),2) - 1;

% KSG Estimation formula.
I = psi(K) + psi(n) - mean(psi(nx+1)+psi(ny+1));
end