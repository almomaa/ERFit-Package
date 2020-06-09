function xdot = lorenzODE(~,x)

xdot = zeros(size(x));

xdot(1) = 10.*(x(2)-x(1));
xdot(2) = 28.*x(1) - x(1).*x(3) - x(2);
xdot(3) = x(1).*x(2) - (8/3).*x(3);
end

