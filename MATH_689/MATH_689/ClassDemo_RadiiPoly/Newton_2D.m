function [x_fixed,y_fixed] = Newton_2D(v0)

lambda = 3;

n=1;N = 5;

x_now = v0(1);
y_now = v0(2);

tol = 1E-15; %Convergence tolerance

while (n <= N+1)
    n
    DF = [8*x_now,  1;
          1,  2*y_now];

    F = [4*x_now^2 + y_now-lambda;
         x_now+y_now^2-1];

    e_now = norm(F, inf)

    h = -DF\F;

    v = [x_now; y_now] + h;

    if (e_now <= 1E-15)
    x_fixed = x_now;
    y_fixed = y_now;
    break;
    end
    
    n=n+1;

    x_now = v(1);
    y_now = v(2);
    %  
end

F = [4*x_now^2 + y_now-lambda;
         x_now+y_now^2-1];


error = norm(F, inf)

%x_now
%y_now
