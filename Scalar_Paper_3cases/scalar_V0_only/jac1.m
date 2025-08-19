 %4D system
function DF = jac1(Y, c, mu, eps)

x = Y(1);
y = Y(2);

%Computing Jacobian
%Uncoment below to get symbolic comps going
%  syms x y c mu eps
% jacobian([x*(1-x)^2;-2*sin(y)*cos(y)*(1-x)+x*((mu+c*(exp(-x.^2)+(eps^2)*exp(-(eps*x).^2)))*(cos(y))^2-(sin(y))^2)], ...
%          [x;y])

DF = [                                                                                            x*(2*x - 2) + (x - 1)^2,                                                                                                            0;
2*cos(y)*sin(y) - sin(y)^2 + cos(y)^2*(mu + (c*eps^2)/cosh(eps*x)) - (c*eps^3*x*sinh(eps*x)*cos(y)^2)/cosh(eps*x)^2, 2*cos(y)^2*(x - 1) - 2*sin(y)^2*(x - 1) - x*(2*cos(y)*sin(y) + 2*cos(y)*sin(y)*(mu + (c*eps^2)/cosh(eps*x)))];
 
end