%
function [y] = fpoly(x,pp,thq);

y = (1/thq(1));
w = y;
for j=2:pp
    w = ((thq(j-1) - x)*w)/thq(j);
    y = y + w;                              
end
