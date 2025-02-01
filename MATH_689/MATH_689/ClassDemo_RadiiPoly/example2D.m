format long
format compact
clear vars

%We are looking at f(x) = (4x_1^2 + x_2 - \lambda; x_1 + x_2^2 - 1)
%Using the Radii Polynomial to validate


lambda = 3;

good = [-0.65;1.28];
%good = [0.79;0.44];
%good = [0.9;-0.3];
%good = [-1.05;-1.43];

[x0,y0] = Newton_2D(good);

p0 = [x0;y0];

[Y0, Z2, p] = radii_polynomial(p0, lambda);

 % Find the roots of p(r) (intersection with y = 0)
 r_roots = roots([Z2, -1, Y0]);
 rminus = min(r_roots)
 rplus = max(r_roots)
 len = rplus;

 intervals = [rminus, rplus, rminus, rplus];

 x_min = p0(1) - len;
 x_max = p0(1) + len;
 y_min = p0(2) - len;
 y_max = p0(2) + len;

 figure(1);
 fplot(p, [-2, 2], 'b-', 'LineWidth', 2);
 hold on;
    
 % Plot the intersection points
 plot(r_roots, zeros(size(r_roots)), 'r*', 'MarkerSize', 8, 'LineWidth', 2);
    
 % Add labels, title, and legend
 xlabel('r');
 ylabel('p(r)');
 title('Radii Polynomial p(r) and Intersection with y = 0');
 grid on;
 legend('p(r)', 'Intersection Points');

 figure(2);
 hold on
 plot(p0(1), p0(2), 'k*','MarkerSize',5)
 rectangle('Position', [x_min, y_min, x_max - x_min, y_max - y_min], ...
              'EdgeColor', 'b', 'LineWidth', 1.5);
% Set the axis limits and labels
xlim([-3, 2]);
ylim([-2.5, 2]);
xlabel('x_1');
ylabel('x_2');
set(gcf,'color','w');
ax = gca;
ax.FontSize = 20; 

hold off;