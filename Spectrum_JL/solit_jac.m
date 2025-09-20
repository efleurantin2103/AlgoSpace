%% Analytical Jacobian of the soliton problem
%
function J = solit_jac(c,ep,x,y)
J = zeros(2,2);
J(1,2) = 1;
J(2,1) = (1-Vep_sol(c,ep,0))-3*(y(1)^2) + Vep_sol(c,ep,x);
end