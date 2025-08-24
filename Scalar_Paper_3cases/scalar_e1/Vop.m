function dydt = Vop(t, Y, options, flag, c, mu, eps, r, u, h)
 
x = Y(1);
y = Y(2);

gamma = mu/(eps^2);

delta = 1e-3*h; % 10*eps; %regularization term

dydt = [x*(1-x)^2 + delta;
         -sin(2*y)*(1-x)+x*((((soliton(r,u,x/(eps*(1-x))))/(eps^2))+gamma+V(c,x,eps))*(cos(y))^2-(sin(y))^2)];
end