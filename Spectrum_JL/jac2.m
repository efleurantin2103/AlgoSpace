function DF = jac2(Y, c, lambda, ep)

x = Y(1);
y = Y(2);

DF = [0,                                      0;
      0, - 2*cos(y)*sin(y) - 2*lambda*cos(y)*sin(y)];
end
 