%% Potential V with magnitude constant c and dependent variable x
%
function f = V1(c,x,eps)
f = c*(eps^2)*sech(eps*(x/(1-x)));
end

