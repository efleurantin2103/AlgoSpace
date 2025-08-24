 %4D system
function DF = jac1(Y, c, mu, eps)

x = Y(1);
y = Y(2);

%Computing Jacobian

% syms x y c mu eps
% jacobian([x*(1-x)^2;-2*sin(y)*cos(y)*(1-x)+x*((mu+c*sech(x))*(cos(y))^2-(sin(y))^2)], ...
%          [x;y])

DF = [                                                                  x*(2*x - 2) + (x - 1)^2,                                                                                                0;
   2*cos(y)*sin(y) - sin(y)^2 + cos(y)^2*(mu + c/cosh(x)) - (c*x*cos(y)^2*sinh(x))/cosh(x)^2, 2*cos(y)^2*(x - 1) - x*(2*cos(y)*sin(y) + 2*cos(y)*sin(y)*(mu + c/cosh(x))) - 2*sin(y)^2*(x - 1)];
 
end