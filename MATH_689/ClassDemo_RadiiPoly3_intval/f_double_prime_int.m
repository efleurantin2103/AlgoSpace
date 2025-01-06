function f_double_prime = f_double_prime_int(x, mu)
    f_double_prime = (intval('6.0')*mu) / ((x - intval('1.0') + mu)^4) + (intval('6.0')*(mu - intval('1.0'))) / ((x + mu)^4);
end