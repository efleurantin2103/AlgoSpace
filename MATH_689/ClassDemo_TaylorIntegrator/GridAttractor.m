%Plotting grid of initial Point evolved to Attractor
clear
format long

tic
numSteps = 10000;
tspan = linspace(0, 150, numSteps); %Time axis

epsilon = 0.25; %Parameters
alpha = 0.95;
gamma = 0.6;
delta = 3.5;
beta = 0.7;
zeta = 0.1;

N = 50;
y0 = [1;1;1];

%Integration:
options=odeset('RelTol',1e-12,'AbsTol',1e-8);
[t,y] = ode45('LangfordField', tspan, y0, options, flag, alpha, beta, delta, gamma, zeta, epsilon);

figure
hold on
    for k = 1:5000
        x0 = y(k,1);
        y0 = y(k,2);
        z0 = y(k,3);
        orbitSeg = LangfordSegment(x0, y0, z0, N, alpha, beta, delta, gamma, zeta, epsilon, 0.04);
        plot3(orbitSeg(1, :), orbitSeg(2, :), orbitSeg(3, :), 'b')
        plot3(x0, y0, z0, 'k*')
    end
 
toc