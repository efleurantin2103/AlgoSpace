%% Potential V with magnitude constant c and dependet variable x
%
function f = V(c,x,eps)
f = c*sech(eps*(x/(1-x)));
end


