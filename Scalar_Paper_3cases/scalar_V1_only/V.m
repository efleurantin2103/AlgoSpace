%% Potential V with magnitude constant c and dependet variable x
% Refer to the manuscript for to comment/uncomment cases as needed
function f = V(c,x,eps)

%f = c*(exp(-(x/(1-x))^2));
%f = c*sech(x/(1-x));
f = c/(1+(x/(1-x))^2)^2;
end

%%%


