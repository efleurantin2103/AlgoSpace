%% Dynamical system formulation of Plucker problem (without singular term)
% for the bvp solver bvp4c
function dydx = bvp_plucker_sys(c,ep,lambda,r,u,x,y)

dydx = zeros(6,1);
dydx(1) = x*(y(2) - y(3)); 
dydx(2) = x*Wm(r,u,x)*y(1) + x*Vep(c,ep,x)*y(1) + x*y(5);
dydx(3) = -x*Wp(r,u,x)*y(1) - x*Vep(c,ep,x)*y(1) - x*y(5); 
dydx(4) = x*lambda*y(1);
dydx(5) = x*Wp(r,u,x)*y(2) + x*Vep(c,ep,x)*y(2) - x*Wm(r,u,x)*y(3) ...
          - x*Vep(c,ep,x)*y(3) - x*2*lambda*y(4);
dydx(6) = x*(1-x)^2;%added this extra line in the compactified system
end