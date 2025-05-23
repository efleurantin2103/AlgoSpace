%Function Aizawa Segment
function orbitSeg = LangfordSegment(x0, y0, z0, N, alpha, beta, delta, gamma, zeta, epsilon, nu)

    a = zeros(1, N+1);
    b = zeros(1, N+1);
    c = zeros(1, N+1);

    a(1) = x0;
    b(1) = y0;
    c(1) = z0;
    a(2) = c(1)*a(1) - beta*a(1) - delta*b(1);
    b(2) = delta*a(1) + c(1)*b(1) - beta*b(1);
    c(2) = gamma + alpha*c(1) - (c(1)*c(1)*c(1))/3 - a(1)*a(1) - b(1)*b(1) ...
        - epsilon*a(1)*a(1)*c(1) - epsilon*b(1)*b(1)*c(1) + zeta*a(1)*a(1)*a(1)*c(1);

for n = 1:N-1
    sum1 = 0;
    sum2 = 0;
    sum3 = 0;
    sum4 = 0;
    sum5 = 0;
    sum6 = 0;
    sum7 = 0;
    sum8 = 0;
    for k = 0:n
        for j = 0:k
            for i = 0:j
                sum1 = sum1 + c(n-k+1)*a(k+1);
                sum2 = sum2 + c(n-k+1)*b(k+1);
                sum3 = sum3 + a(n-k+1)*a(k+1);
                sum4 = sum4 + b(n-k+1)*b(k+1);
                sum5 = sum5 + c(n-k+1)*c(k-j+1)*c(j+1);
                sum6 = sum6 + a(n-k+1)*a(k-j+1)*c(j+1);
                sum7 = sum7 + b(n-k+1)*b(k-j+1)*c(j+1);
                sum8 = sum8 + a(n-k+1)*a(k-j+1)*a(j-i+1)*c(i+1);
            end
        end
    end

a(n+2) = (sum1 - beta*a(n+1) - delta*b(n+1))/(n+1);
b(n+2) = (delta*a(n+1) + sum2 - beta*b(n+1))/(n+1);
c(n+2) = (alpha*c(n+1) - (sum5)/3 - sum3 - sum4 - epsilon*sum6 - epsilon*sum7 + zeta*sum8)/(n+1);
end

numPoints = 1000;
T = linspace(-nu, nu, numPoints);
orbitPoints = zeros(3, numPoints);

    for j = 1:numPoints
        this_t = T(j);
        thisPoint = [0;0;0];
         for n = 0:N
            thisCoef = [a(n+1); b(n+1); c(n+1)];
            thisPoint = thisPoint + thisCoef*this_t^n;
         end
        orbitPoints(1:3, j) = thisPoint;
    end

orbitSeg = orbitPoints;