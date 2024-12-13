function DF = jac3(Y, c, mu, eps)

x = Y(1);
y = Y(2);

%%% uncomment to compute Jacobian for a given potential

% syms x y c mu eps
% jacobian([x*(1-x)^2;
%          -2*sin(y)*cos(y)*(1-x)+x*((mu+c*(eps^2)*exp(-(eps*(x)).^2))*(cos(y))^2-(sin(y))^2)], ...
%          [x;y])

DF = [                                                                                      x*(2*x - 2) + (x - 1)^2,                                                                                                         0;
     2*cos(y)*sin(y) - sin(y)^2 + cos(y)^2*(mu + c*eps^2*exp(-eps^2*x^2)) - 2*c*eps^4*x^2*exp(-eps^2*x^2)*cos(y)^2, 2*cos(y)^2*(x - 1) - x*(2*cos(y)*sin(y) + 2*cos(y)*sin(y)*(mu + c*eps^2*exp(-eps^2*x^2))) - 2*sin(y)^2*(x - 1)];
end
