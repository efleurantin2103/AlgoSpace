%%%%% Function Differential.m
function DF =Differential(x, mu)

% % Define the symbolic variables
% syms x mu
% 
% % Define the function f(x)
% f_sym = x + (mu / ((x - 1 + mu)^2)) + ((mu - 1) / ((x + mu)^2));
% 
% % Compute the derivative of f(x) with respect to x
% df_sym = diff(f_sym, x)

DF = 1 - (2*(mu - 1))/(mu + x)^3 - (2*mu)/(mu + x - 1)^3;

%%%%%%