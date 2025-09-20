%% Boundary conditions for soliton
%
function bc = solit_bc(xmax,ya,yb)
bc = [ya(2);
       (yb(1) + (xmax/(1 + xmax))*yb(2))];
end