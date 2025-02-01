function [Y0, Z0, Z1, Z2, p] = radii_polynomial(p,nu)

a0=p(1);

Dfp0 = Differential(p);
fp0 = F(p);

A = inv(Dfp0);

Y0 = calc_Y0(p,A*fp0,nu,2);

Z0 = compute_norm(eye(3)-A*Dfp0, nu, 2);

Z1 = calc_Z1(p, nu, 2);

% Calculate the terms inside the max function
term1 = compute_norm(A, nu, 2);
term2 = 1/(2*abs(a0));
    
% Calculate Z2 using the max function
Z2 = 2*max(term1, term2);
    
 % Define the resulting radii polynomial p(r) 
 p = @(r) Z2 * r.^2 - (1 - Z0 - Z1)*r + Y0;

end