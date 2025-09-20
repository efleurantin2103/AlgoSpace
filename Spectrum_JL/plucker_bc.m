%% Initial condition for pluckers system (coded as a boundary condition)
%
function bc = plucker_bc(ya,yb)
bc = zeros(6,1);
bc(1,1) = ya(1)-1;
bc(2,1) = ya(2);
bc(3,1) = ya(3);
bc(4,1) = ya(4);
bc(5,1) = ya(5);
bc(5,1) = ya(6);
end