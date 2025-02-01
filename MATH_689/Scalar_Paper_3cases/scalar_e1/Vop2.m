function out = Vop2(t, Y, options, flag, c, mu, eps, r, u)
 
x = Y(1);
y = Y(2);

% Or use more robust form:
safex = min(max(x, 0), 0.99999);

out = [x*(1-x)^2;
         -2*sin(y)*cos(y)*(1-x)+x*((soliton(r,u,safex/(1-safex)) + mu +V1(c,safex,eps))*(cos(y))^2-(sin(y))^2)];
