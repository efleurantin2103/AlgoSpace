function DF =Differential_intval(x, mu)

% % Define the symbolic variables
% syms x mu
% 
% % Define the function f(x)
% f_sym = x + (mu / ((x - 1 + mu)^2)) + ((mu - 1) / ((x + mu)^2));
% 
% % Compute the derivative of f(x) with respect to x
% df_sym = diff(f_sym, x)

DF = intval('1.0') - (intval('2.0')*(mu - 1))/(mu + x)^3 - (intval('2.0')*mu)/(mu + x - intval('1.0'))^3;
