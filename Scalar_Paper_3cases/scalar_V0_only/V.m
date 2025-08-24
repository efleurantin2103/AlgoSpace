%% Potential V with magnitude constant c and dependet variable x
% refer to manuscript and uncomment/comment for different cases
function f = V(c,x,eps)

%f = c*(exp(-(x/(1-x))^2));
f = c/((1+2*(x/(1-x))^4));
%f = c/((cosh(1.2*(x/(1-x))))^2);
end

%%%
