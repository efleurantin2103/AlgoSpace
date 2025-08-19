
function DF = jac2(Y, c, mu, eps)

x = Y(1);
y = Y(2);

%Computing Jacobian

% 
% syms x y c mu
% jacobian([x*(1-x)^2;-2*sin(y)*cos(y)*(1-x)+x*((mu+(c/(1+x^2)))*(cos(y))^2-(sin(y))^2)], ...
%          [x;y])

DF = [                                                                  x*(2*x - 2) + (x - 1)^2,                                                                                                  0;
2*cos(y)*sin(y) - sin(y)^2 + cos(y)^2*(mu + c/(x^2 + 1)) - (2*c*x^2*cos(y)^2)/(x^2 + 1)^2, 2*cos(y)^2*(x - 1) - 2*sin(y)^2*(x - 1) - x*(2*cos(y)*sin(y) + 2*cos(y)*sin(y)*(mu + c/(x^2 + 1)))];
 
end