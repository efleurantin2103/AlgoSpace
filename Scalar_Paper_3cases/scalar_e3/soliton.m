%% Interpolated soliton function
%
function f = soliton(x,y,z)
s=max(x);
if (z > s)
    f = 0;
else
    f = interp1(x,y,z);
end
