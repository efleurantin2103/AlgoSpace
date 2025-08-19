%Code for our ODE
function dydt = Vop(t, Y, options, flag, c, mu, eps)
 
x = Y(1);
y = Y(2);

dydt = [ x*(1-x)^2;
         -2*sin(y)*cos(y)*(1-x)+x*((mu+V(c,x,eps))*(cos(y))^2-(sin(y))^2)];