function a = accessFourierCoeff(n, A, N)

if (-N <=n) && (n <= N)
    a  = A(N+1+n);
else
    a = 0;
end

return