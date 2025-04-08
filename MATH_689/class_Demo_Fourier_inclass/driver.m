clear all
format long

N = 50;
A = emptyFourier(N);
B = A;

g = @(x) (1/(1+sin(x)^2));

for n = -N:N
    %n
    cn = computeComplexFourierCoefficients(g,n);
    B = setFourierCoef(n, cn, B, N);
    %A = setFourierCoef(n, 1/(n^4+1), A, N);
end

numPoints = 1000;
x_axis = linspace(0, 2*pi, numPoints);
f_values = zeros(1, numPoints);
g_values = zeros(1, numPoints);
for k = 1: numPoints
    x = x_axis(k);
    f_values(k) = evaluateFourierSeries(x,B,N);
    g_values(k) = g(x);
end

figure
hold on
plot(x_axis, f_values, 'r')
plot(x_axis, g_values, 'g')

x = 1
error = abs(evaluateFourierSeries(x,B,N)-g(x))

%As far as HW goes, you may use or modify this code
%You may need to also code a function
%multiplyFourierCoef(A, B, N)
