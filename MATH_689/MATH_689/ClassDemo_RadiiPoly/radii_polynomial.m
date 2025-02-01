function [Y0, Z2, p] = radii_polynomial(x, lambda)

x1=x(1);
x2=x(2);

Dfp0 = Differential(x, lambda);
fp0 = F(x, lambda);

A = inv(Dfp0);

Y0 = norm(A*fp0, inf);

    % Calculate the terms inside the max function
    term1 = (16*abs(x2) + 2) / (abs(16*x1*x2 - 1));
    term2 = (8 + 16*abs(x1)) / (abs(16*x1*x2 - 1));
    
    % Calculate Z2 using the max function
    Z2 = max(term1, term2);
    
    % Define the resulting radii polynomial p(r) 
    p = @(r) Z2 * r.^2 - r + Y0;
end