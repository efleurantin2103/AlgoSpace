function F =F(p)

a0 = p(1);
a1 = p(2);
a2 = p(3);

F = [a0^2 - (1/3);...
     2*a0*a1 - 1;...
    a2*a0+a1*a1+a0*a2];
