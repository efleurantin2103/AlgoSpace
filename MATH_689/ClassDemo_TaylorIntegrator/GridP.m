clear
format long

epsilon = 0.25; %Parameters
alpha = 0.95;
gamma = 0.6;
delta = 3.5;
beta = 0.7;
zeta = 0.1;

N = 50;

numPoints = 7;

the_xs = linspace(-2, 2, numPoints);
the_ys = linspace(-2, 2, numPoints);
the_zs = linspace(-2, 2, numPoints);

figure
hold on
for i = 1:numPoints
    i
    numPoints
        for j = 1:numPoints
            for k = 1:numPoints
                x0 = the_xs(i);
                y0 = the_ys(j);
                z0 = the_zs(k);
                orbitSeg = LangfordSegment(x0, y0, z0, N, alpha, beta, delta, gamma, zeta, epsilon, 0.04);
                plot3(orbitSeg(1, :), orbitSeg(2, :), orbitSeg(3, :), 'b')
                plot3(x0, y0, z0, 'k*')
            end
        end
end