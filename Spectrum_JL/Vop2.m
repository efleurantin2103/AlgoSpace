%Code for our ODE
function dydt = Vop2(t, Y, options, flag, c, lambda, ep)
 
x = Y(1);
y = Y(2);

gamma = (1-lambda)/(ep^2);


dydt = [ x*(1-x)^2;
         -sin(2*y)*(1-x)+x*((gamma+c*exp(-(x/(1-x))^2))*(cos(y))^2-(sin(y))^2)];