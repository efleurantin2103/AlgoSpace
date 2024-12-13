clear all
format long

%Langford Parameters:
alpha = 0.95; beta = 0.7; delta = 3.5; gamma = 0.6; zeta = 0.1; epsilon = 0.25;

%Equilibria:
p0 = [0; 0; -1.1050]; 
p1 = [0; 0; -0.83824]; 
p2 = [0; 0; 1.9433];

%Choose an equilibrium point to compute manifolds:
p = p1;

Dfp = langfordDifferential(p, ...
      epsilon, alpha, gamma, delta, beta, zeta);

[Q, Lambda] = eigs(Dfp);

%The names of the eigenvalues have to be
%set by hand:
lambda_u = Lambda(3,3);
xi_u = Q(:, 3);

lambda_s1 = Lambda(1,1);
lambda_s2 = Lambda(2,2);
xi_s1 = Q(:, 1);
xi_s2 = Q(:, 2);

N = 35;

scale = 0.75;%Be careful of this one

P = zeros(N+1, N+1, 3);

%Compute a two dimensional manifold:
SM = computeStableM(p,Dfp,lambda_s1, lambda_s2, xi_s1, xi_s2, N, scale,epsilon, alpha, gamma, delta, beta, zeta);

numR = 20;
numTheta = 50;

the_rs = linspace(0.001, 1, numR);
the_theta = linspace(0, 2*pi, numTheta);
numPoints = numR*numTheta;

Ws_loc = zeros(3, numPoints);

pointNum = 1;
for m = 1:numR
    m
    numR
      for n = 1:numTheta
        this_r = the_rs(m);
        this_theta = the_theta(n);
        sigma1 = this_r*cos(this_theta);
        sigma2 = this_r*sin(this_theta);
        thisPoint = [0; 0; 0];
            for order = 0:N
                for vader = 0:order
                    m1 = order-vader;
                    m2 = vader;

                    thisPoint = thisPoint + ...
                    reshape(SM(m1+1, m2+1, :), [3,1])*...
                    ((sigma1 + i*sigma2)^m1)*((sigma1 - i*sigma2)^m2);
                end
            end
        Ws_loc(:, pointNum) = real(thisPoint);
        pointNum = pointNum + 1;
    end
end

figure
hold on
plot3(Ws_loc(1, :)', Ws_loc(2, :)', Ws_loc(3, :)', 'b.')
plot3(p0(1), p0(2), p0(3), 'b*')
plot3(p1(1), p1(2), p1(3), 'r*')
plot3(p2(1), p2(2), p2(3), 'g*')

%%%%%%%%%%%%% Plotting global manifold %%%%%%%%%%%%%%%%
numPoints = 50;
theTheta = linspace(0, 2*pi, numPoints);
Ws = zeros(3, 500*numPoints);
count = 1;
for k = 1:numPoints
    k

    thisTheta = theTheta(k);
    the_x = cos(thisTheta);
    the_y = sin(thisTheta);

    bestPoint = evaluate_complexCase(the_x, the_y, ...
    reshape(SM(:, :, 1), [N+1, N+1]), ...
    reshape(SM(:, :, 2), [N+1, N+1]), ...
    reshape(SM(:, :, 3), [N+1, N+1]), N);

    T = 1.6;

    tspan = linspace(0, -T, 500);
    IC = bestPoint';

    %Integration:
    options=odeset('RelTol',1e-13,'AbsTol',1e-13);
    [t, guessOrbit] = ode45('langfordField', tspan, ...
    IC, options, flag, ...
    alpha, beta, delta, gamma, zeta, epsilon);

    Ws(:, 500*(k-1)+1:500*k) = guessOrbit';
    
    count = count + 1;
end

figure
hold on
plot3(Ws_loc(1, :)', Ws_loc(2, :)', Ws_loc(3, :)', 'b.')
plot3(Ws(1, :)', Ws(2, :)', Ws(3, :)', 'k.')
plot3(p0(1), p0(2), p0(3), 'b*')
plot3(p1(1), p1(2), p1(3), 'r*')
