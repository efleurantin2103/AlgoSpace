function F =F_intval(x, mu)

F = x + (mu / ((x - intval('1.0') + mu)^2)) + ((mu - intval('1.0')) / ((x + mu)^2));
