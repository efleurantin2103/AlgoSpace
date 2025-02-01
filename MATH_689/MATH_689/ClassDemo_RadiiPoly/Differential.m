function DF =Differential(p, lambda)

x = p(1);
y = p(2);

DF = [8*x,  1;
      1,  2*y];
