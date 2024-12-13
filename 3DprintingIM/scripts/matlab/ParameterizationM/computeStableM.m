%Function Compute Stable Manifold
function SM = computeStableM(p,Dfp,lambda_s1, lambda_s2, xi_s1, xi_s2, ...
                        N, scale,epsilon, alpha, gamma, delta, beta, zeta)

P = zeros(N+1, N+1, 3);

P(1, 1, :) = p;
P(2, 1, :) = scale*xi_s1;
P(1, 2, :) = scale*xi_s2;

Id = eye(3);

coeffDecay = zeros(1, N+1);
coeffDecay(1) = norm(P(:, 1), inf);
coeffDecay(2) = max([norm(reshape(P(1, 2, :), [3, 1]), inf), ...
norm(reshape(P(2, 1, :),[3, 1]), inf)]);

for order = 2:N
  order
  theseAbs = [];
   for steps = 0:order
    n1 = order - steps;
    n2 = steps;
    A = Dfp - (n1*lambda_s1 + n2*lambda_s2)*Id;
    sum1 = 0;
    sum2 = 0;
    sum3 = 0;
    sum4 = 0;
    sum5 = 0;
    sum6 = 0;
    sum7 = 0;
    sum8 = 0;
        for k1 = 0:n1
            for k2 = 0:n2
                sum1 = sum1 + P(n1-k1+1, n2-k2+1, 3)* P(k1+1, k2+1, 1);
                sum2 = sum2 + P(n1-k1+1, n2-k2+1, 3)* P(k1+1, k2+1, 2);
                sum3 = sum3 + P(n1-k1+1, n2-k2+1, 1)*P(k1+1, k2+1, 1);
                sum4 = sum4 + P(n1-k1+1, n2-k2+1, 2)*P(k1+1, k2+1, 2);
                 for j1 = 0:k1
                     for j2 = 0:k2
                        sum5 = sum5 + P(n1-k1+1, n2-k2+1, 3)* P(k1-j1+1, k2-j2+1, 3)* P(j1+1, j2+1, 3);
                        sum6 = sum6 + P(n1-k1+1, n2-k2+1, 1)* P(k1-j1+1, k2-j2+1, 1)* P(j1+1, j2+1, 3);
                        sum7 = sum7 + P(n1-k1+1, n2-k2+1, 2)* P(k1-j1+1, k2-j2+1, 2)* P(j1+1, j2+1, 3);
                         for i1 = 0:j1
                            for i2 = 0:j2
                                sum8 = sum8 + P(n1-k1+1, n2-k2+1, 1)* P(k1-j1+1, k2-j2+1, 1)* P(j1-i1+1, j2-i2+1, 1)*P(i1+1, i2+1, 3);
                            end
                         end
                     end
                 end
             end
         end
    sum1 = sum1 - P(n1+1, n2+1, 3)*P(1, 1, 1) - P(n1+1, n2+1, 1)*P(1, 1, 3) ;
    sum2 = sum2 - P(n1+1, n2+1, 3)*P(1, 1, 2) - P(n1+1, n2+1, 2)*P(1, 1, 3);
    sum3 = sum3 - P(n1+1, n2+1, 1)*P(1, 1, 1) - P(n1+1, n2+1, 1)*P(1, 1, 1);
    sum4 = sum4 - P(n1+1, n2+1, 2)*P(1, 1, 2) - P(n1+1, n2+1, 2)*P(1, 1, 2);
    sum5 = sum5 - 3*P(1,1,3)^2*P(n1+1, n2+1, 3);
    sum6 = sum6 - 2*P(1,1,1)*P(1,1,3)*P(n1+1, n2+1, 1) - (P(1,1,1)^2)*P(n1+1, n2+1, 3);
    sum7 = sum7 - 2*P(1,1,2)*P(1,1,3)*P(n1+1, n2+1, 2) - (P(1,1,2)^2)*P(n1+1, n2+1, 3);
    sum8 = sum8 - 3*P(1,1,1)^2*P(1,1,3)*P(n1-j1+1, n2-j2+1, 1) - ...
    P(1,1,1)^3*P(n1-j1+1, n2-j2+1, 3);

    Rmn = [ -sum1;
    -sum2;
    (sum5)/3 + sum3 + sum4 + ...
    epsilon*sum6 + epsilon*sum7 ...
    - zeta*sum8];

    P(n1+1, n2+1, :) = A\Rmn;

    theseAbs = [theseAbs, ...
    norm(reshape(P(n1+1, n2+1, :), [3 1]), inf)];
   end
  coeffDecay(order+1) = max(theseAbs);
end

SM = P;

figure
hold on
plot([0:N], log(coeffDecay)/log(10))
