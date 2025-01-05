%%%%%%% Function F.m
function F =F(x, mu)

F = x + (mu / ((x - 1 + mu)^2)) + ((mu - 1) / ((x + mu)^2));
%%%%%%%%
