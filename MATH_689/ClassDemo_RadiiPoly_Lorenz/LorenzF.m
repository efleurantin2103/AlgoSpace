function F = lorenzF(Y, sigma, rho, beta)

x = Y(1);
y = Y(2);
z = Y(3);


F = [sigma*(y-x); 
    x*(rho - z) - y;
    x*y - beta*z];
