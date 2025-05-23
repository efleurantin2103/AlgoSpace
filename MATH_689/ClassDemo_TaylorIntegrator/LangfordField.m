%Function aizawa field
function dydt = LangfordField(t, Y, options, flag, alpha, beta, delta, gamma, zeta, epsilon)

x = Y(1);
y = Y(2);
z = Y(3);

dydt = [(z-beta)*x - delta*y;
        delta*x + (z - beta)*y;
        gamma + alpha*z - (z*z*z)/3 - (x*x + y*y)*(1 + epsilon*z) + zeta*z*x*x*x];