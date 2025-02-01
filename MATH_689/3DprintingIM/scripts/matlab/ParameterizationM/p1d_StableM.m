clear all
format long

alpha = 0.95; beta = 0.7; delta = 3.5; gamma = 0.6; zeta = 0.1; epsilon = 0.25;

N = 50;

a = zeros(1, N+1);
b = zeros(1, N+1);
c = zeros(1, N+1);

a(1) = 0;
b(1) = 0;
c(1) = 1.9433;

scale = 3*10^(-2);

a(2) = scale*0;
b(2) = scale*0;
c(2) = scale*1;
p0 = [a(1);b(1);c(1)];
p1 = [a(2);b(2);c(2)];

for n = 2:N
    I = [1 0 0;
        0 1 0;
        0 0 1];

    DF = [1.2433 -3.5 0;
          3.5 1.2433 0;
          0 0 -2.09254];

    E = -2.09254;

    A = DF-n*E*I;

    B = inv(A);

    sum1 = 0;
    sum2 = 0;
    sum3 = 0;
    sum4 = 0;
    sum5 = 0;
    sum6 = 0;
    sum7 = 0;
    sum8 = 0;
    for k = 1:n-1
        sum1 = sum1 + c(n-k+1)*a(k+1);
        sum2 = sum2 + c(n-k+1)*b(k+1);
        sum3 = sum3 + a(n-k+1)*a(k+1);
        sum4 = sum4 + b(n-k+1)*b(k+1);
        for j = 0:k
            sum5 = sum5 + c(n-k+1)*c(k-j+1)*c(j+1);
            sum6 = sum6 + a(n-k+1)*a(k-j+1)*c(j+1);
            sum7 = sum7 + b(n-k+1)*b(k-j+1)*c(j+1);
            for i = 0:j
                sum8 = sum8 + a(n-k+1)*a(k-j+1)*a(j-i+1)*c(i+1);
            end
        end
    end
    a(n+1) = - sum1;
    b(n+1) = - sum2;
    c(n+1) = (sum5)/3 + sum3 + sum4 + ...
    epsilon*sum6 + epsilon*sum7 - ...
    zeta*sum8;

    S = [a(2:n+1); b(2:n+1); c(2:n+1)];

    S1 = B*S;

    S2 = transpose(S1);

    for x = 2:n
    m{x} = S1(:,x);
    end
end
steps = 40;

numPoints = 100;

thetas = linspace(-1, 1, numPoints);

m0 = transpose(p0);
m1 = transpose(p1);

P =[m0; m1; S2];

Ws = zeros(3, numPoints);

for k = 1:numPoints
    thisTheta = thetas(k);
    value = [0;0;0];
    for j = 0:N
        value = value + P(j+1, :)'*thisTheta^j;
    end
    Ws(:, k) = value;
end

%ROTATE PICTURE IN 3D
figure
hold on
plot3(Ws(1, :), Ws(2, :), Ws(3, :), 'b')
plot3(a(1), b(1), c(1), 'k*')
