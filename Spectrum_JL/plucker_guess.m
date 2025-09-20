%% Initial Guess function for Plucker system
%
function f = plucker_guess(x)
f = zeros(6,1);
f(1,1) = 1/(sqrt(1+x^2));
f(2,1) = -(x/(1+x^2)^1.5);
f(3,1) = -(x/(1+x^2)^1.5);
f(4,1) = -(x/(1+x^2)^1.5);
f(5,1) = -(x/(1+x^2)^1.5);
f(6,1) = -(x/(1+x^2)^1.5); %added this extra line here for the 
                           %compatified system. Not sure if this the right
                           %one
end
