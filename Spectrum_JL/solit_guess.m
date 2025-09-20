%% Initial guess function for soliton
%
function f = solit_guess(amp0,x)
f = [(amp0*sech(x))/sqrt(1 + x^2);
     (-amp0*((sech(x)*(x + (1 + x^2)*tanh(x)))/(1 + x^2)^1.5))];
end
