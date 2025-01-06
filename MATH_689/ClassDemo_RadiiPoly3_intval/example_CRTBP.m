%Dependencies: The Codes Require the IntLab library in order to run

format long
format compact
clear vars

%Using the Radii Polynomial to validate

%%%%%%%%%% add new path

  wng = warning;          % supress warning that "format" is overloaded
  warning off
  addpath('Intlab_V13')

  path(path)			         % make sure paths are correct
								 % works in all Matlab versions                       
  warning(wng);
  startup   %start Intlab

%%%%%%%%%%

good = Newton_1D(-0.83);

mu = intval('1.0')/(intval('1.0123'));

p0 =  intval(good);

[Y0, Z0, Z2, p] = radii_polynomial(p0, mu);

%return 
% Find the roots of p(r) (intersection with y = 0)
r_roots = roots([sup(Z2), -(1-sup(Z0)), sup(Y0)]);
rminus = min(r_roots)
rplus = max(r_roots)
len = rplus;

 new_p = @(r) sup(Z2) * r.^2 - (1-sup(Z0))*r + sup(Y0);
 x_min = good - len;
 x_max = good + len;

 figure(1);
 fplot(new_p, [-0.05, 0.1], 'b-', 'LineWidth', 2);
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

figure(2)
hold on
plot(good, 0, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
% Plot the interval line
plot([x_min, x_max], [0, 0], 'k-', 'LineWidth', 2);
% Plot the brackets
text(x_min, 0, '[', 'FontSize', 24, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
text(x_max, 0, ']', 'FontSize', 24, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
xlim([x_min-0.01, x_max+0.01]);
%ylim([-0.01, 0.01]);
xlabel('Interval');
title('Interval of Existence');
% Remove y-axis tick labels
set(gca, 'YTick', []);
% Display the plot
grid on;
hold off