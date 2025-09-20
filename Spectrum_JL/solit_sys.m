%% Dynamical systems formulation for computing soliton (without singular term)
%
function f = solit_sys(c,ep,x,y)
f = [y(2);
        (1-Vep_sol(c,ep,0))*y(1) - y(1)^3 + Vep_sol(c,ep,x)*y(1)];
end
