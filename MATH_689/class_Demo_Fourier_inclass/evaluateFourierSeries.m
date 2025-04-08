function fx = evaluateFourierSeries(x, A, N)

sum = 0;

for n = -N:N
    sum = sum + accessFourierCoeff(n, A, N)*exp(1i*n*x);
end

if abs(imag(sum)) > 10^(-6)
    crashNow = crashHere
end

fx = real(sum);

return