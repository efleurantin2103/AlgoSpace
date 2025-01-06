function x_fixed = Newton_1D(v0)

n=1;N = 6;

mu = 1/(1.0123);

x_now = v0;

tol = 1E-15; %Convergence tolerance

while (n <= N+1)
    n
    DF = Differential(x_now, mu);

    F = x_now + (mu / ((x_now - 1 + mu)^2)) + ((mu - 1) / ((x_now + mu)^2));

    e_now = norm(F, inf)

    h = -DF\F;

    v = x_now + h;

    if (e_now <= 1E-15)
    x_fixed = x_now;
    break;
    end
    
    n=n+1;

    x_now = v;
    %  
end

F = x_now + (mu / ((x_now - 1 + mu)^2)) + ((mu - 1) / ((x_now + mu)^2))


error = norm(F, inf)

%x_now
%y_now
