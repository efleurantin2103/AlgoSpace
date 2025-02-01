function F =F(p, lambda)

x = p(1);
y = p(2);

F = [4*x^2 + y-lambda;
         x+y^2-1];
