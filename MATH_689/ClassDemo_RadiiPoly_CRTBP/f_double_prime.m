%%%%%% Function f_double_prime.m
function f_double_prime = f_double_prime(x, mu)
    f_double_prime = (6*mu) / ((x - 1 + mu)^4) + (6*(mu - 1)) / ((x + mu)^4);
end
%%%%%%