%Function Evaluate Complex Case
function f=evaluate_complexCase(v, w, A, B, C, N)

z1 = (v + i*w);
z2 = (v - i*w);

Ax = zeros(N+1);
By = Ax;
Cz = Ax;

for m = 1:(N+1)
    thisCoef = A(1:(N+1-(m-1)), m);
    Ax(m) = horner(z1, thisCoef, N-(m-1));

    thisCoef = B(1:(N+1-(m-1)), m);
    By(m) = horner(z1, thisCoef, N-(m-1));
    
    thisCoef = C(1:(N+1-(m-1)), m);
    Cz(m) = horner(z1, thisCoef, N-(m-1));
end

Poly_x = horner(z2, Ax, N);
Poly_y = horner(z2, By, N);
Poly_z = horner(z2, Cz, N);

f=real([ Poly_x Poly_y Poly_z]);