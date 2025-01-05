function DF=lorenzDifferential(p, sigma, rho, beta)

x = p(1);
y = p(2);
z = p(3);

DF=[  -sigma,      sigma, 0;
     rho - z, -1,        -x;
     y,        x,     -beta];