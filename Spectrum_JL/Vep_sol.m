%% Epsilon scaled potential Vep created from potential V
%
function f = Vep_sol(c,ep,x)
f = (ep^2)*V_sol(c,ep*x);
end