%% Entropic Regression
%
% 

%% Citation
%  To cite this code:
%  [ 1 ] Abd AlRahman R. AlMomani, Jie Sun, and Erik Bollt. How Entropic 
%        Regression Beats the Outliers Problem in Nonlinear System 
%        Identification. Chaos 30, 013107 (2020).
%        
%% Description
% S = erfit(Phi, Xdot): solve the inverse problem:
% given Xdot and Phi, find S such that Xdot = Phi*S.
% erfit returns a matrix S of coefficient for a linear regression of 
% the responses in vector field Xdot as a linear combination of the 
% basis functions in the matrix Phi. 
% 


%% Inputs
%%
%
% * Phi  : Nxd basis matrix (i.e. polynomial expansion of the state variable).
% * Xdot : Nxr vector field (derivative of the state variables).
% * options: Structure that has the erfit options. See eroptset.m
%

%% Outputs
% S : dxr coefficients matrix.
% Info: structure that contain different infomation about the process.

%% 
% 
% Author: Abd AlRahman R. AlMomani
% Clarkson University, 2020.
% Version 1.1.0
% For any comments and feedback please email the author at:
% aaalmoma@clarkson.edu, or, almomaniar@gmail.com

%% Function Body
%
%%
function [S, Info] = erfit(Phi, Xdot, varargin)

err = 'All imputs should have the same number of rows.';
assert((size(Phi,1)==size(Xdot,1)), err);

if (~isempty(varargin)) && (isstruct(varargin{end}))
    options = varargin{end};   
    
elseif exist('eroptset','file')
    options = eroptset;
end

dim = size(Xdot,2); [m,N] = size(Phi);


if options.useparallel
    pp = gcp;
    options.NumWorkers = pp.NumWorkers;
else
    options.NumWorkers = 0;
end


%Initialize the solution
S = zeros(N,dim);


% Start Entropic Regression
for i=1:dim
    % default tolerence initialization
    options.tol = 0;
    if ~options.EmbeddedShuffleTest
        In = zeros(options.numPerm,1);
        for j=1:options.numPerm
            In(j) = options.MIEstimator(Xdot(:,i),Xdot(randperm(m),i),options);
        end
        options.tol = quantile(In,options.alpha);
    end

    Info.estimatedTolerence(i) = options.tol;

    
    % For the ith dimension:
    IX = 1:N; %Explore all the candidate functions
    success = false;
    if ~options.skipForward
        [IX1, success] = erforward(Phi,Xdot(:,i),IX,options);
    end

    if (~success || options.skipForward)
        IX1 = IX; 
    end

    Info.ForwardSelectedIndex{i} = IX1;
    Info.ForwardSuccess(i) = success;
    
    

    % Eliminate the weak candidate functions from IX
    % through the backward ER.
    IX2 = erbackward(Phi(:,IX1), Xdot(:,i), 1:length(IX1),options);
    
    
    IX = IX1(IX2);
    Info.BackwardEliminatedIndex{i} = setdiff(IX1,IX);
    
    % re-check for the constant term (add it to the estimated indices)
    if ~ismember(1,IX), IX = cat(2,1,IX); end
    
    % Estimate parameters magnitude
    index = false(size(S,1),1); index(IX) = true;
    S(:,i) = options.grayModelEstimator(Phi,Xdot(:,i),index);
    
    if ((abs(S(1,i))<options.h) && (sum(~~S(:,i))>1))
        S(1,i) = 0;
        IX = setdiff(IX,1);
    end
    
    Info.Index{i} = IX;
    Info.S{i} = S(IX,i);
end

end



