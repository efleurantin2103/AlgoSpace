%% Potential V with magnitude constant c and dependet variable x
%
function f = V(c,x,eps)
%f = c*(x^2)*exp(-x);
%change x to x/(1-x)
%f = c*(exp(-(x/(1-x)).^2)+(eps^2)*exp(-(eps*x/(1-x)).^2));
%f = c*(eps^2)*exp(-(eps*x/(1-x)).^2)+1;
%f = c*(eps^2)*exp(-eps^2*(x/(1-x))^2)+eps^2;
f = c*exp(-((x/(1-x))).^2);
end

%%%
%Slow fast scales and see that the winding occurs seperately one on the
%fast scale and one on the slow scale
