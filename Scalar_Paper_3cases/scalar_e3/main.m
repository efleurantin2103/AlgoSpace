clear all;
format long;

eps=0.1;
c=-30;

%eigenvalue paramaters
mu=eps^2*linspace(0.01,10,10);

%Compute V_0
r = linspace(0,1000,238); %Create r vector from 0 to 1000 
u = -3*((1)./(cosh(1.2*r))).^2;% Calculate u


lim=1;
scale2=10^(-5);

for h = 1:1 %length(mu)

%Getting fixed points
gamma = ((mu(h))/eps^2);
Fpb=atan(-sqrt(gamma));
Fpa=atan(sqrt(gamma));
Fpbb=atan(-sqrt(mu(h)));
Fpaa=atan(sqrt(mu(h)));
theta=-pi;
theta2=pi;
theta3=2*pi;
eq1 = [0;0];
eq2 = [0;theta];
eq3 = [0;theta2];
eq4 = [0;theta3];

Fp1=[1;Fpb-2*pi];
Fp2=[1;Fpa-2*pi];
Fpp1=[1;Fpbb-pi];
Fpp2=[1;Fpaa-pi];

%Jacobians
DFl = jac3(eq1, c, mu(h), eps);
DF = jac3(eq2, c, mu(h), eps);
DF3 = jac3(eq3, c, mu(h), eps);
DF4 = jac3(eq4, c, mu(h), eps);

[Ql, Lambdal] = eig(DFl);
[Q, Lambda] = eig(DF);
[Q3, Lambda3] = eig(DF3);
[Q4, Lambda4] = eig(DF4);

%eigenvalues
lambda_1l = Lambdal(1,1);
lambda_2l = Lambdal(2,2);
lambda_1 = Lambda(1,1);
lambda_2 = Lambda(2,2);
lambda_13 = Lambda3(1,1);
lambda_23 = Lambda3(2,2);
lambda_14 = Lambda4(1,1);
lambda_24 = Lambda4(2,2);

%eigenvectors
xi_1l = Ql(:, 1);
xi_2l = Ql(:, 2);
xi_1 = Q(:, 1);
xi_2 = Q(:, 2);
xi_13 = Q3(:, 1);
xi_23 = Q3(:, 2);
xi_14 = Q4(:, 1);
xi_24 = Q4(:, 2);

DF1 = jac4(Fp2, c, mu(h), eps);

% Compute the eigenvectors and eigenvalues of the Jacobian for Wuc
[V1, D] = eig(DF1);

diag(D)
V1(:,1)
V1(:,2)

% Extract the center/stable eigenvector
if D(1,1) < D(2,2)    
    V1 = V1(:,2);
else
    V1 = V1(:,1);
end

scale = 10^(-4);

x0 = Fp2 - scale*V1;
x1 = eq2 + scale*xi_2;
x3 = eq3 + scale*xi_23;
x4 = eq4 + scale*xi_24;
xl = eq1 + scale*xi_2l;

%Integrate
options=odeset('RelTol',1e-13,'AbsTol',1e-13);
options2=odeset('RelTol',(1e-8)*h^(-5/2),'AbsTol',1e-13);
[t, W0u2l] = ode45('Vop2',[0 10000], xl, options, flag, c, mu(h), eps, r, u);
[t, W0u2] = ode45('Vop2',[0 10000], x1, options, flag, c, mu(h), eps, r, u);
[t, W0u23] = ode45('Vop2',[0 10000], x3, options, flag, c, mu(h), eps, r, u);
[t, W0u24] = ode45('Vop2',[0 10000], x4, options, flag, c, mu(h), eps, r, u);
%[t, W0c] = ode45('Vop',[0 -1E9], x0, options, flag, c, mu(h), eps, r, u);
%[t, W0c] = ode45('Vop',[2e4 0], x0, options, flag, c, mu(h), eps, r, u, h);
[t, W0c] = ode15s('Vop',[2e4 0], x0, options2, flag, c, mu(h), eps, r, u, h);



if h == 1
    line_color = 'k';  % Black for first index
    line_width = 2.5;  % Thicker for black
else
    line_color = 'r';  % Red for other indices
    line_width = 0.5;  % Thicker for black
end

figure(2)

% Left subplot
subplot(1, 2, 1)
hold on
plot(W0u2l(:, 1), W0u2l(:, 2), [line_color '-'], 'LineWidth', line_width)
%plot(W0u2(:, 1), W0u2(:, 2), 'r-', 'LineWidth', 0.5)
plot(W0u23(:, 1), W0u23(:, 2), 'r-', 'LineWidth', 0.5)
plot(W0u24(:, 1), W0u24(:, 2), 'r-', 'LineWidth', 0.5)
plot(Fpp1(1),Fpp1(2)-2*pi,'b.', 'MarkerSize',20)
plot(Fpp1(1),Fpp1(2)+pi,'b.', 'MarkerSize',20)
plot(Fpp1(1),Fpp1(2)-pi,'b.', 'MarkerSize',15)
plot(Fpp1(1),Fpp1(2),'b.', 'MarkerSize',15)
plot(Fpp1(1),Fpp1(2)+2*pi,'b.', 'MarkerSize',15)
plot(Fpp2(1),Fpp2(2)-pi,'r.', 'MarkerSize',20)
plot(Fpp2(1),Fpp2(2)-2*pi,'r.', 'MarkerSize',20)
plot(Fpp2(1),Fpp2(2)+pi,'r.', 'MarkerSize',15)
plot(Fpp2(1),Fpp2(2),'r.', 'MarkerSize',15)
plot(Fpp2(1),Fpp2(2)+2*pi,'r.', 'MarkerSize',15)
plot(eq1(1),eq1(2),'r.', 'MarkerSize',15)
%plot(eq2(1),eq2(2),'r.', 'MarkerSize',15)
plot(eq3(1),eq3(2),'r.', 'MarkerSize',15)
plot(eq4(1),eq4(2),'r.', 'MarkerSize',15)
xlabel('$\sigma$', 'Interpreter', 'latex')
ylabel('$\theta$', 'Interpreter', 'latex')
yyaxis left
axis([0 1 -10 9])
yyaxis right
axis([0 1 -10 9])
set(gca, 'YColor', 'k','FontSize', 20, 'TickLabelInterpreter', 'latex','YLabel', [])
hold off

% Right subplot
subplot(1, 2, 2)
hold on
plot(W0c(:, 1), W0c(:, 2), [line_color '-'], ...
     'LineWidth', line_width, ...
     'LineStyle', '-', ...
     'LineJoin', 'round', ...
     'Marker', 'none')
plot(Fp2(1),Fp2(2),'b.', 'MarkerSize',15)
plot(W0c(:, 1), -3*pi/2 + 0*W0c(:, 1),W0c(:, 1), -pi/2 + 0*W0c(:, 1),W0c(:, 1), pi/2 + 0*W0c(:, 1),W0c(:, 1), 3*pi/2 + 0*W0c(:, 1))
% plot(Fp2(1),Fp2(2)+pi,'b.', 'MarkerSize',15)
% plot(Fp2(1),Fp2(2)+2*pi,'b.', 'MarkerSize',15)
% plot(Fp1(1),Fp1(2),'r.', 'MarkerSize',15)
%plot(Fp1(1),Fp1(2)+pi,'r.', 'MarkerSize',15)
% plot(Fp1(1),Fp1(2)+2*pi,'r.', 'MarkerSize',15)
plot(0,3*pi/2,'b.', 'MarkerSize',15)
plot(0,pi/2,'b.', 'MarkerSize',15)
plot(0,-pi/2,'b.', 'MarkerSize',15)
plot(0,5*pi/2,'b.', 'MarkerSize',15)
yyaxis left
axis([0 1 -10 9])
set(gca, 'YColor', 'k','YLabel', [])
yyaxis right
axis([0 1 -10 9])
set(gca, 'YColor', 'k', 'FontSize', 20, 'TickLabelInterpreter', 'latex')
xlabel('$\tau$', 'Interpreter', 'latex')
ylabel('$\psi$', 'Interpreter', 'latex')
hold off
floor((W0c(end,2)-W0c(1,2))/pi)


% With additional formatting
sgtitle('Scenario 3: $V_{1,\varepsilon}(x)=-\frac{30\varepsilon^2}{(1+(\varepsilon x)^2)^2}$,  $V_0(x) = -\frac{3}{\cosh^2(1.2x)}$', 'Interpreter', 'latex', ...
    'FontSize', 20, ...
    'FontWeight', 'bold')

% Final touches
set(gcf, 'Color', 'w')
drawnow

end
<<<<<<< Updated upstream
=======

print(2, '-depsc', '-painters', 'scalar3_case1')
>>>>>>> Stashed changes
