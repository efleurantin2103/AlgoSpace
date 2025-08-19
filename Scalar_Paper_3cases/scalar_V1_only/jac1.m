 %4D system
function DF = jac1(Y, c, mu, eps)

x = Y(1);
y = Y(2);

%Computing Jacobian

%  syms x y c mu eps
% jacobian([x*(1-x)^2;-2*sin(y)*cos(y)*(1-x)+x*((mu+c*(exp(-x.^2)))*(cos(y))^2-(sin(y))^2)], ...
%          [x;y])

% DF = [                                                                                                                                                                                                        x*(2*x - 2) + (x - 1)^2,                                                                                                                       0;
%       2*cos(y)*sin(y) - sin(y)^2 + cos(y)^2*(mu + (c*x^2*exp(x/(x - 1)))/(x - 1)^2) - x*cos(y)^2*((2*c*x^2*exp(x/(x - 1)))/(x - 1)^3 - (2*c*x*exp(x/(x - 1)))/(x - 1)^2 + (c*x^2*exp(x/(x - 1))*(x/(x - 1)^2 - 1/(x - 1)))/(x - 1)^2), 2*cos(y)^2*(x - 1) - 2*sin(y)^2*(x - 1) - x*(2*cos(y)*sin(y) + 2*cos(y)*sin(y)*(mu + (c*x^2*exp(x/(x - 1)))/(x - 1)^2))];


% syms x y c mu eps
% jacobian([x*(1-x)^2;
%          -2*sin(y)*cos(y)*(1-x)+x*((mu+c*eps^2*sech(eps*x))*(cos(y))^2-(sin(y))^2)], ...
%          [x;y])

% DF = [                                                                                x*(2*x - 2) + (x - 1)^2,                                                                                                    0;
%       2*cos(y)*sin(y) + cos(y)^2*(mu + c*x^2*exp(-x)) - sin(y)^2 + x*cos(y)^2*(2*c*x*exp(-x) - c*x^2*exp(-x)), 2*cos(y)^2*(x - 1) - x*(2*cos(y)*sin(y) + 2*cos(y)*sin(y)*(mu + c*x^2*exp(-x))) - 2*sin(y)^2*(x - 1)];

DF = [                                                              x*(2*x - 2) + (x - 1)^2,                                                                                                  0;
2*cos(y)*sin(y) - sin(y)^2 + cos(y)^2*(mu + c*exp(-x^2)) - 2*c*x^2*exp(-x^2)*cos(y)^2, 2*cos(y)^2*(x - 1) - 2*sin(y)^2*(x - 1) - x*(2*cos(y)*sin(y) + 2*cos(y)*sin(y)*(mu + c*exp(-x^2)))];
 
end