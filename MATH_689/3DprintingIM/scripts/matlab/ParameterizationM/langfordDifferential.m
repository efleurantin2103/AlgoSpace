%Function Langford Differential
%1) To compute a jacobian
%2) To get the linear part of the homological equations
function DF = langfordDifferential(p, ...
                epsilon, alpha, gamma, delta, beta, zeta)

DF = [p(3)-beta, -delta, p(1);
delta, p(3)-beta, p(2);
-2*p(1)*(1+epsilon*p(3))+3*zeta*p(3)*p(1)*p(1), ...
-2*p(2)*(1+epsilon*p(3)), ...
alpha-p(3)*p(3)-epsilon*(p(1)*p(1)+p(2)*p(2))+zeta*p(1)*p(1)*p(1)];
end