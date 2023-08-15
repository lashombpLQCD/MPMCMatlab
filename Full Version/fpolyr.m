%
function [y] = fpoly(x,p,thq);

xp = 1;
for i=1:p
    xp = xp - (1/thq(i))*x*xp;
end
y = 1 - xp;