
function DF = jac2(Y, c, mu, eps)

x = Y(1);
y = Y(2);

%Computing Jacobian
%Uncomment below to get symbolic comps going
% syms x y c mu
% jacobian([x*(1-x)^2;-2*sin(y)*cos(y)*(1-x)+x*((mu+(c/((1+2*x)^4)))*(cos(y))^2-(sin(y))^2)], ...
%          [x;y])

DF = [                                                                  x*(2*x - 2) + (x - 1)^2,                                                                                                    0;
2*cos(y)*sin(y) - sin(y)^2 + cos(y)^2*(mu + c/(2*x + 1)^4) - (8*c*x*cos(y)^2)/(2*x + 1)^5, 2*cos(y)^2*(x - 1) - 2*sin(y)^2*(x - 1) - x*(2*cos(y)*sin(y) + 2*cos(y)*sin(y)*(mu + c/(2*x + 1)^4))];
 
end