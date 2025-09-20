%% Analytical Jacobian of boundary conditions for soliton problem
%
function [dBCdya, dBCdyb] = solit_bc_jac(xmax,ya,yb)
    dBCdya = [0 1; 0 0];
    dBCdyb = [0 0; 1 (xmax/(1 + xmax))];
end