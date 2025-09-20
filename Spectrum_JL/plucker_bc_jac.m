%% Analytical Jacobian of boundary conditions for Plucker problem 
%
function [dBCdya, dBCdyb] = plucker_bc_jac(ya,yb)
    dBCdya = eye(6);
    dBCdyb = zeros(6,6);
end