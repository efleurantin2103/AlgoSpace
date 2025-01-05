%%%%%% Function radii_polynomial.m
function [Y0, Z0, Z2, p] = radii_polynomial(x, mu)

Dfp0 = Differential(x, mu);
fp0 = F(x, mu);

A = inv(Dfp0);

Y0 = norm(A*fp0, inf);

Z0 = norm(1 - A*Dfp0, inf);

r_star = 0.01;

Z2 = A * max(abs(f_double_prime(x-r_star, mu)));
    
% Define the resulting radii polynomial p(r) 
 p = @(r) Z2 * r.^2 - (1-Z0)*r + Y0;
end