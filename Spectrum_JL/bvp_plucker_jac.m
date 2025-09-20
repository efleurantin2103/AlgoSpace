%% Analytic Jacobian of the Dynamical system bvp_plucker_sys.m
%
function J = bvp_plucker_jac(c,ep,lambda,r,u,x,y)
J = zeros(6,6);
J(1,2) = x;
J(1,3) = -x;
J(1,6) = y(2)-y(3);
J(2,1) = x*(Wm(r,u,x) + Vep(c,ep,x));
J(2,5) = x;
J(2,6) = y(1)*(Wm(r,u,x) + Vep(c,ep,x))+y(5);
J(3,1) = x*(-Wp(r,u,x)- Vep(c,ep,x));
J(3,5) = -x;
j(3,6) = y(1)*(-Wp(r,u,x)- Vep(c,ep,x))-y(5);
J(4,1) = x*lambda;
J(4,6) = y(1)*lambda;
J(5,2) = x*(Wp(r,u,x)+ Vep(c,ep,x));
J(5,3) = x*(-Wm(r,u,x) - Vep(c,ep,x));
J(5,4) = -2*lambda*x;
J(5,6) = y(2)*(Wp(r,u,x)+ Vep(c,ep,x)) + y(3)*(-Wm(r,u,x) - Vep(c,ep,x)) -2*lambda*y(4);
J(6,6) = (1-x)^2 - 2*x*(1-x);
end