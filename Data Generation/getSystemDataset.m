function [varargout] = getSystemDataset(odefun,varargin)

varargout = cell(1,nargout);

p = inputParser;

validateFunHanle  = @(x) isa(x,'function_handle');
validPosScalar    = @(x) isscalar(x) && isnumeric(x) && (x>0);
validPosScalarLT1 = @(x) isscalar(x) && isnumeric(x) && (x>=0) && (x<=1);
addRequired(p,'odefun',validateFunHanle);

p.addParameter('odeSolver', @ode45,validateFunHanle);
p.addParameter('SampleSize', 1000, validPosScalar);

p.addParameter('derivativeEstimator',  @centralDifference);

p.addParameter('dim', 3, validPosScalar);

p.addParameter('tao', 0.01, validPosScalar);

p.addParameter('eps1', 0.005, @isscalar);

p.addParameter('eps2', 0, @isscalar);

p.addParameter('Berp', 0, validPosScalarLT1);

p.addParameter('expanOrder', 4, validPosScalar);

validfun = @(x) (isnumeric(x) || strcmp(x,'rand'));
p.addParameter('initCondition', 'rand',validfun);


p.addParameter('skipTrans', true);


parse(p,odefun,varargin{:});
options = p.Results;


odefun = options.odefun;

if isnumeric(options.initCondition)
    x0 = options.initCondition;
else
    x0 = rand(options.dim,1);
end            
Info.x0 = x0;

h = options.tao;
N = options.SampleSize;

X = [];
if options.skipTrans
    [~,X] = options.odeSolver(odefun, [0 100], x0); %Transient
end
if ~isempty(X)
    x0 = X(end,:)';
end
Info.x0 = x0;
[t,X] = options.odeSolver(odefun, 0:h:(N+2)*h, x0); %Data Sampling

Info.Xclean = X;
Info.t      = t;


X = X + options.eps1*randn(size(X));

Info.Xbasenoise = X;


IX = (rand(size(X,1),1)<=options.Berp);
if sum(IX)
    X(IX,:) = X(IX,:) + options.eps2*randn(size(X(IX,:)));
end

Info.Xcorrupted = X;
Info.CorruptionIndex = IX;

[Xdot, X] = options.derivativeEstimator(X,h);

Info.Xdot = Xdot;

[Phi, P] = polyspace(X,options.expanOrder);


Info.Phi = Phi;
Info.PowrMatrix = P;
Info.h = h;
Info.eps1 = options.eps1;
Info.eps2 = options.eps2;
Info.Berp = options.Berp;
Info.ExpansionOrder = options.expanOrder;

if nargout==1
    varargout{1} = Info;
elseif nargout==2
    varargout{1} = Phi;
    varargout{2} = Xdot;
elseif nargout==3
    varargout{1} = Phi;
    varargout{2} = Xdot;
    varargout{3} = Info;
end

end

