clearvars; close all; clc;

x = [2.50582594582818,4.79803200521955,8.04276503609942];
% x = rand(1,3);

P = [-10,-13.3326397311082,0;
      10, 10.1272728412525,0.986596914100953;
      -0.1,-0.5,-1];

for i=1:5000
    plot3(x(1),x(2),x(3),'*b')
    hold on
    drawnow
    
    xd = x*P;
    x1 = 0.01*xd + x.*(0.01*randi(1,3)+0.99) ;
    x = x1;
%     x = x1 - 0.5*pdist2(x,x1);
    
end