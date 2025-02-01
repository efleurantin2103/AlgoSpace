format long
format compact
clear vars

%Using the Radii Polynomial to validate
 nu=1/4;

 p0=[0.57735026918962;0.86602540378443;-0.64951905283832];

 [Y0, Z0, Z1, Z2, p] = radii_polynomial(p0,nu);

% Find the roots of p(r) (intersection with y = 0)
 r_roots = roots([Z2, - (1 - Z0 - Z1), Y0]);
 rminus = min(r_roots)
 rplus = max(r_roots)

 figure(1);
 fplot(p, [-0.15, 0.4], 'b-', 'LineWidth', 2);
 hold on;
    
 % Plot the intersection points
 plot(r_roots, zeros(size(r_roots)), 'r*', 'MarkerSize', 8, 'LineWidth', 2);

 % Add labels, title, and legend
 xlabel('r');
 ylabel('p(r)');
 title('Radii Polynomial p(r) and Intersection with y = 0');
 grid on;
 legend('p(r)', 'Intersection Points');
 hold off

% getting the function

% Set the range and number of points
lam_min = 0.1;
lam_max = 0.6;

points = 1000;

% Generate points along the u-axis
lam_vals = linspace(lam_min, lam_max, points);

% Solve the equation for v
x_N = p0(3)*(lam_vals-(1/3)).^2+p0(2)*(lam_vals-(1/3))+p0(1);
x_true = sqrt(lam_vals);

% Create the formatted string
str = sprintf('sup |x^(N)(\\lambda) - x(\\lambda)| \\leq %.14f,\n|', rminus);

% Display the string
disp(str);


figure(2);
hold on
plot(1/3, sqrt(1/3), 'k*','MarkerSize',10)
plot(lam_vals, x_N, 'b-', 'LineWidth', 2);
plot(lam_vals, x_true, 'r-', 'LineWidth', 2);
xlabel('\lambda');
ylabel('x(\lambda)');
set(gcf,'color','w');
axis([0 0.65 -inf inf])
ax = gca;
ax.FontSize = 20; 
hold off;

%Comments:
% Figure(2) is the plot of  \bar{x}^(N)(\lambda) for N = 2 show as a blue curve versus
% and \bar{x}(\lambda) = sqrt{\lambda} shown in red. The point (\lambda_0; x0) = (1/3; \sqrt{1=3}) 
% is the black star. The approximation error is rigorously proven 
% to be less than 0.045 in the indicated domain.
