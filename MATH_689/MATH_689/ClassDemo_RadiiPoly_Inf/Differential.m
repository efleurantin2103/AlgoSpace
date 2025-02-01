function DF =Differential(p)

a0 = p(1);
a1 = p(2);
a2 = p(3);

%%%%%%%%%%%%%%%%%%%%%%%%
%%Uncomment below for matlab to compute symbolic Jacobian
%
% syms a0 a1 a2
% % %
% jacobian([a0^2 -(1/3), ...
%           2*a0*a1 - 1, ...
%           a2*a0+a1*a1+a0*a2],[a0,a1,a2])
%%%%%%%%%%%%%%%%%%%%%%%%%%


DF = [2*a0,    0,    0;
      2*a1, 2*a0,    0;
      2*a2, 2*a1, 2*a0];
