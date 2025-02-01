%% Potential V with magnitude constant c and dependet variable x
%
function f = V(c,x,eps)
f = c*(1/(1+(eps*(x/(1-x)))^2));
end

%%%
%Slow fast scales and see that the winding occurs seperately one on the
%fast scale and one on the slow scale
