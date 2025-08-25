clear all;
format long;
%This code computes all three cases, so one has to be careful and
%comment/uncomment selected cases. Please refer to the manuscript.

eps=0;%Not needed   
c=-2.6;%This is one is for scenario 2 
mu=0;%eigenvalue parameter  

lim=1;%set amount of Fp

for i=1:lim
theta=(i-1)*pi;
eq2 = [0;theta];%Fp at sigma=0

DF = jac1(eq2, c, mu, eps);

[Q, Lambda] = eig(DF);

%assigning values
lambda_1 = Lambda(1,1);
lambda_2 = Lambda(2,2);

xi_1 = Q(:, 1);
xi_2 = Q(:, 2);

scale = 10^(-4);
x1 = eq2 + scale*xi_2;

options=odeset('RelTol',1e-13,'AbsTol',1e-13);
[t, W0u2] = ode45('Vop',[0 5000], x1, options, flag, c, mu, eps);

figure(1)
hold on
plot(eq2(1), eq2(2), 'r.', 'MarkerSize', 15)
plot(W0u2(:, 1), W0u2(:, 2), 'r', 'LineWidth', 3)
xline(1,'k')
set(gca, 'YColor', 'k', 'FontSize', 20, 'TickLabelInterpreter', 'latex')
xlabel('$\sigma$', 'Interpreter', 'latex')
ylabel('$\theta$', 'Interpreter', 'latex')
hold off

set(gcf, 'Color', 'w')

end

%print(1, '-depsc', '-painters', 'scalar_eigen_case3')
% 
