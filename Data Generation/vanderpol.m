function xdot = vanderpol(~,x)

xdot = [x(2);...
        5*(1-x(1).^2).*x(2) - x(1)];
end

