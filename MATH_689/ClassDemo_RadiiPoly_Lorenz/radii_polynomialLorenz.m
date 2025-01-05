function [Y0, Z0, Z2, p] = radii_polynomialLorenz(x, sigma, rho, beta)

x1=x(1);
x2=x(2);
x3=x(3);

Dfp0 = lorenzDifferential(x, sigma, rho, beta);
fp0 = LorenzF(x, sigma, rho, beta);

A = inv(Dfp0);

I = diag(ones(3,1));

Z0 = norm(I-A*Dfp0, inf); 

Y0 = norm(A*fp0, inf);

% Extract the elements of A
a12 = A(1,2);
a13 = A(1,3);
a21 = A(2,1);
a23 = A(2,3);
a31 = A(3,1);
a32 = A(3,2);

% Calculate the terms inside the max function
term1 = 2*(abs(a12) + abs(a13));
term2 = 2*(abs(a21) + abs(a23));
term3 = 2*(abs(a31) + abs(a32));

% Calculate Z2 using the max function
Z2 = max([term1, term2, term3]);
    
   
 % Define the resulting radii polynomial p(r) 
 p = @(r) Z2 * r.^2 - (1 - Z0)*r + Y0;
end