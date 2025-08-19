function DF = jac4(Y, c, mu, eps)

x = Y(1);
y = Y(2);

%% uncomment to compute Jacobian for a given potential

% syms x y c mu eps
% jacobian([x*(1-x)^2;
%          -sin(2*y)*(1-x)+x*(((mu)/eps^2)*(cos(y))^2-(sin(y))^2)], ...
%          [x;y])



DF = [                  x*(2*x - 2) + (x - 1)^2,                                                                0;
sin(2*y) - sin(y)^2 + (mu*cos(y)^2)/eps^2, 2*cos(2*y)*(x - 1) - x*(2*cos(y)*sin(y) + (2*mu*cos(y)*sin(y))/eps^2)];

end
