%%%% example_CRTBP.m

format long
format compact
clear vars

p0 = Newton_1D(-0.83);

mu = 1/(1.0123);

[Y0, Z0, Z2, p] = radii_polynomial(p0, mu);

% Find the roots of p(r) (intersection with y = 0)
r_roots = roots([Z2, -1, Y0]);
rminus = min(r_roots)
rplus = max(r_roots)
len = rplus;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%Conclusion
%if r \in [rminus, rplus] \cap [0, r_star] = [rminus, r_star] then p(r)<0.
%We conclude the existence of a unique \bar{x} \in [p0 - rminus, p0 + rminus]. 
%%%%%%%%%%%%%%%%%%%%%%%%%%

x_min = p0 - len;
x_max = p0 + len;

 figure(1);
 fplot(p, [-0.05, 0.1], 'b-', 'LineWidth', 2);
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
lam_min = -10;
lam_max = 10;

points = 10000;

% Generate points along the u-axis
lam_vals = linspace(lam_min, lam_max, points);

% Solve the equation for v
fx = lam_vals + (mu ./ ((lam_vals - 1 + mu).^2)) + ((mu - 1) ./ ((lam_vals + mu).^2));

figure(2);
hold on
plot(p0, 0, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
plot(lam_vals, fx, 'b-', 'LineWidth', 2);
xlabel('x');
ylabel('f(x)');
set(gcf,'color','w');
axis([-11 11 -11 11])
%xlim([x_min-0.01, x_max+0.01]);
ax = gca;
ax.FontSize = 20; 
% Plot the brackets
text(x_min, 0, '[', 'FontSize', 24, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
text(x_max, 0, ']', 'FontSize', 24, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
% Add x-axis and y-axis
xline(0, 'k--', 'LineWidth', 1.5); % x-axis
yline(0, 'k--', 'LineWidth', 1.5); % y-axis
hold off


%%%%%%%%%%%%