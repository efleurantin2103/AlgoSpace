clear all;
format long;

eps=0.1;
c=-30;

mu=linspace(0.05,0.25,10 ...
    );
%4.347388643264107
%% 
r = linspace(0,1000,238); %Create r vector from 0 to 1000 
u = -(0.1*exp(-r)).^2;% Calculate u

% Create the plot
% figure
% plot(r, u, 'LineWidth', 1.5)
% xlabel('r')
% ylabel('u')
%% 

Points2=[];

lim=1;
scale2=10^(-5);

for h = 1:length(mu)
gamma = ((mu(h))/eps^2);
Fpb=atan(-sqrt(gamma));
Fpa=atan(sqrt(gamma));
Fpbb=atan(-sqrt(mu(h)));
Fpaa=atan(sqrt(mu(h)));
theta=pi;
theta2=2*pi;
theta3=3*pi;
eq1 = [0;0];
eq2 = [0;theta];
eq3 = [0;theta2];
eq4 = [0;theta3];

% Fp1=[1;Fpb-2*pi];
Fp1=[1;Fpb-2*pi];
Fp2=[1;Fpa-2*pi];
Fpp1=[1;Fpbb-pi];
Fpp2=[1;Fpaa-pi];

DFl = jac3(eq1, c, mu(h), eps);
DF = jac3(eq2, c, mu(h), eps);
DF3 = jac3(eq3, c, mu(h), eps);
DF4 = jac3(eq4, c, mu(h), eps);

[Ql, Lambdal] = eig(DFl);
[Q, Lambda] = eig(DF);
[Q3, Lambda3] = eig(DF3);
[Q4, Lambda4] = eig(DF4);

lambda_1l = Lambdal(1,1);
lambda_2l = Lambdal(2,2);
lambda_1 = Lambda(1,1);
lambda_2 = Lambda(2,2);
lambda_13 = Lambda3(1,1);
lambda_23 = Lambda3(2,2);
lambda_14 = Lambda4(1,1);
lambda_24 = Lambda4(2,2);

xi_1l = Ql(:, 1);
xi_2l = Ql(:, 2);
xi_1 = Q(:, 1);
xi_2 = Q(:, 2);
xi_13 = Q3(:, 1);
xi_23 = Q3(:, 2);
xi_14 = Q4(:, 1);
xi_24 = Q4(:, 2);

DF1 = jac4(Fp2, c, mu(h), eps);

% Compute the eigenvectors and eigenvalues of the Jacobian
[V1, D] = eig(DF1);

diag(D)
V1(:,1)
V1(:,2)
%return
% Extract the center/stable eigenvector
if D(1,1) < D(2,2)    
    V1 = V1(:,2);
else
    V1 = V1(:,1);
end

scale = 10^(-4);
x0 = Fp2 - scale*V1;%-scale2;
x1 = eq2 + scale*xi_2;
x3 = eq3 + scale*xi_23;
x4 = eq4 + scale*xi_24;
xl = eq1 + scale*xi_2l;

options=odeset('RelTol',1e-13,'AbsTol',1e-13);
[t, W0u2l] = ode45('Vop2',[0 10000], xl, options, flag, c, mu(h), eps, r, u);
[t, W0u2] = ode45('Vop2',[0 10000], x1, options, flag, c, mu(h), eps, r, u);
[t, W0u23] = ode45('Vop2',[0 10000], x3, options, flag, c, mu(h), eps, r, u);
[t, W0u24] = ode45('Vop2',[0 10000], x4, options, flag, c, mu(h), eps, r, u);
[t, W0c] = ode45('Vop',[0 -1E7], x0, options, flag, c, mu(h), eps, r, u);

alpha=0.75;
threshold = eps^(alpha);
threshold2 = 1-eps^(alpha);%This works better and matches the neq parameterization in Section 5
cross_indices = find(diff(W0c(:,1) >= threshold) ~= 0);
cross_indices2 = find(diff(W0u2(:,1) >= threshold2) ~= 0);
cross_indices3 = find(diff(W0u23(:,1) >= threshold2) ~= 0);
cross_indices4 = find(diff(W0u24(:,1) >= threshold2) ~= 0);
cross_indicesl = find(diff(W0u2l(:,1) >= threshold2) ~= 0);

theta_diff = abs(W0c(cross_indices,2) - W0c(1,2));
theta_diff2 = abs(W0u2(cross_indices2,2) - W0u2(1,2));
theta_diff3 = abs(W0u23(cross_indices3,2) - W0u23(1,2));
theta_diff4 = abs(W0u24(cross_indices4,2) - W0u23(1,2));

figure(2)
set(gcf, 'Color', 'w', 'Position', [100, 100, 600, 400])

% Left subplot
subplot(1, 2, 1)
hold on
plot(W0u2l(1:cross_indicesl, 1), W0u2l(1:cross_indicesl, 2), 'r-', 'LineWidth', 0.5)
plot(W0u2(1:cross_indices2, 1), W0u2(1:cross_indices2, 2), 'r-', 'LineWidth', 0.5)
plot(W0u23(1:cross_indices3, 1), W0u23(1:cross_indices3, 2), 'r-', 'LineWidth', 0.5)
plot(eq1(1),eq1(2),'r.', 'MarkerSize',15)
plot(eq2(1),eq2(2),'r.', 'MarkerSize',15)
plot(eq3(1),eq3(2),'r.', 'MarkerSize',15)
xline(threshold2,'--', 'LineWidth', 1.5);
xlabel('$\sigma$', 'Interpreter', 'latex')
ylabel('$\theta$', 'Interpreter', 'latex')
yyaxis left
axis([0 1 -8 8])
yyaxis right
axis([0 1 -8 8])
set(gca, 'YColor', 'none','FontSize', 20, 'TickLabelInterpreter', 'latex','YLabel', [],'YTick', [])
hold off

% Right subplot
subplot(1, 2, 2)
hold on
plot(W0c(1:cross_indices, 1), W0c(1:cross_indices, 2), 'r', ...
     'LineWidth', .5, ...
     'LineStyle', '-', ...
     'LineJoin', 'round', ...
     'Marker', 'none')
xline(threshold,'--', 'LineWidth', 1.5);
plot(Fp2(1),Fp2(2),'r.', 'MarkerSize',15)
yyaxis left
axis([0 1 -8 8])
set(gca, 'YColor', 'none', 'YTick', [])
yyaxis right
axis([0 1 -8 8])
set(gca, 'YColor', 'k', 'FontSize', 20, 'TickLabelInterpreter', 'latex')
xlabel('$\tau$', 'Interpreter', 'latex')
ylabel('$\psi$', 'Interpreter', 'latex')
hold off

% Final touches
set(gcf, 'Color', 'w')
drawnow

end
