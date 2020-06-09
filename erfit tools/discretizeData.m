function [Symbols,Edges] = discretizeData(X,options)

% Discretize Data
%
%   This function is part of the Entropic Regression package (erfit)
%   
%
% Inputs:
%    X: is an n-by-d matrix. n points in the d-dimensional sample space.
%    options: Options structure that have the fields:
%      'BinMethod': string that have the name of the method used to 
%                   discretize the data. The default values for
%                   'BinMethod' is 'fixed', which will use fixed number 
%                   of bins that is determined by the option:
%      'numBins': The default value for 'numBins' is 16.
%
%      For all possible BinMethod that can be used, see the options
%      function ( eroptset ).
 
% Outputs:
%    Symbols: is an n-by-d integer matrix. (Labels).
%    Edges: The edges values of the intervals that is used to label the
%    data.
%
%
% ---------------------------------------------------
%   Abd AlRahman R. AlMomani
%   aaalmoma(at)clarkson(dot)edu
%   May 2020
% ---------------------------------------------------

if strcmp(options.BinMethod,'fixed')
    [~,Edges] = histcounts(X,options.numBins);
else
    [~,Edges] = histcounts(X,'BinMethod',options.BinMethod);
end
Symbols  = discretize(X,Edges);
end

