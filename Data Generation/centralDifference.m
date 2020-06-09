function [Xdot, X] = centralDifference(X,h)
% Derivative using central difference method
Xdot = (1/(2*h))*(X(3:end,:)-X(1:end-2,:));
X(1,:) = []; X(end,:) = [];
end

