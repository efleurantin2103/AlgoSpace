% Set up the figure with professional styling
set(0, 'defaultAxesFontSize', 12);
set(0, 'defaultTextFontSize', 12);
set(0, 'defaultAxesFontName', 'Times New Roman');
set(0, 'defaultTextFontName', 'Times New Roman');

% Define the system parameters
x = 0:0.01:1;
y = -1.5:0.01:1;
[X, Y] = meshgrid(x, y);

% Define the vector field
Ydot = -Y.*(1+X) - (1-X).*Y.^2;
Xdot = X.*(X-1);

% Create figure with specific size
fig = figure('Position', [100, 100, 800, 400], 'Color', 'white');

% Set axes background to white
%set(gca, 'Color', 'white');

% Create streamslice plot with customized appearance
h = streamslice(X, Y, Xdot, Ydot, 0.15);
set(h, 'Color', [0.2 0.2 0.2], 'LineWidth', 0.8);

% Customize grid and axis appearance
grid on;
set(gca, 'GridAlpha', 0.15, 'LineWidth', 1);
box on;
axis equal;
axis([-0.15 1 -1.5 1]);

% Add fixed points with improved markers
hold on;
plot(0, -1, 'o', 'Color', [0.8 0 0], 'MarkerSize', 7, 'LineWidth', 3);
plot(0, 0, '.', 'Color', [0.8 0 0], 'MarkerSize', 30, 'LineWidth', 1.5);

% Add the β̄=0 line
plot([-0.15 1], [0 0], '--', 'Color', [0 0.5 0], 'LineWidth', 2);

% Add the s=0 line
plot([0 0], [-1.5 1], '--', 'Color', [0 0.5 0], 'LineWidth', 2);

% Label axes with LaTeX formatting
xlabel('$s$', 'Interpreter', 'latex', 'FontSize', 30);
ylabel('$\bar{\beta}$', 'Interpreter', 'latex', 'FontSize', 26);

% Add title
%title('Phase Portrait', 'FontSize', 14);

% Add legend
%legend('', 'Stable Node', 'Unstable Node', '$\bar{\beta}=0$', 'Interpreter', 'latex', 'Location', 'best');

% Adjust overall plot appearance
set(gca, 'TickLength', [0.02 0.02]);
set(gca, 'Layer', 'top');

% Adjust spacing
set(gcf, 'PaperPositionMode', 'auto');
tight = get(gca, 'TightInset');
set(gca, 'Position', [tight(1)+0.05 tight(2)+0.05 0.95-tight(1)-tight(3) 0.95-tight(2)-tight(4)]);

% Save the figure in high resolution
%print(fig, 'endp_flow', '-dpdf', '-r300');
%print(fig, 'endp_flow', '-dsvg', '-r300');
print(fig, 'endp_flow', '-depsc', '-painters', '-r300');