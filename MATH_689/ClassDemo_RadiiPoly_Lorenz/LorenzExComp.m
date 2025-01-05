format long
format compact
clear vars

%Lorenz Parameters
sigma = 10;
beta = 8/3;
rho = 28;   

%Equilibrium Solns:
%p0 = [0;0;0];
p0 = [-sqrt(beta*(rho-1));-sqrt(beta*(rho-1));rho-1];
%p0 = [-sqrt(beta*(rho-1));-sqrt(beta*(rho-1));rho-1];

[Y0, Z0, Z2, p] = radii_polynomialLorenz(p0, sigma, rho, beta);

 % Find the roots of p(r) (intersection with y = 0)
 r_roots = roots([Z2, -(1 - Z0), Y0]);
 rminus = min(r_roots)
 rplus = max(r_roots)
 

 figure(1);
 fplot(p, [-2, 6], 'b-', 'LineWidth', 2);
 hold on;
    
 % Plot the intersection points
 plot(r_roots, zeros(size(r_roots)), 'r*', 'MarkerSize', 8, 'LineWidth', 2);
    
 % Add labels, title, and legend
 xlabel('r');
 ylabel('p(r)');
 title('Radii Polynomial p(r) and Intersection with y = 0');
 grid on;
 legend('p(r)', 'Intersection Points');

 %len = abs(rplus - rminus)/2;
 len = rplus;

 x_min = p0(1) - len;
 x_max = p0(1) + len;
 y_min = p0(2) - len;
 y_max = p0(2) + len;
 z_min = p0(3) - len;
 z_max = p0(3) + len;

 % Create the vertices of the box
vertices = [x_min, y_min, z_min;
            x_max, y_min, z_min;
            x_max, y_max, z_min;
            x_min, y_max, z_min;
            x_min, y_min, z_max;
            x_max, y_min, z_max;
            x_max, y_max, z_max;
            x_min, y_max, z_max];

% Define the edges of the box
edges = [1, 2; 2, 3; 3, 4; 4, 1; 1, 5; 2, 6; 3, 7; 4, 8; 5, 6; 6, 7; 7, 8; 8, 5];


 figure(2);
 hold on
 plot3(p0(1), p0(2), p0(3),'k*','MarkerSize',5)
% Plot the edges of the box
for i = 1:size(edges, 1)
    start_vertex = vertices(edges(i, 1), :);
    end_vertex = vertices(edges(i, 2), :);
    plot3([start_vertex(1), end_vertex(1)], [start_vertex(2), end_vertex(2)], [start_vertex(3), end_vertex(3)], 'b-', 'LineWidth', 1.5);
end

% Set the axis limits and labels
xlim([p0(1)-15, p0(1)+15]);
ylim([p0(2)-15, p0(2)+15]);
zlim([p0(3)-15, p0(3)+15]);
xlabel('X');
ylabel('Y');
zlabel('Z');

% Set the aspect ratio to make the box look proportional
%daspect([1, 1, 1]);

% Display the plot
hold off;
grid on;