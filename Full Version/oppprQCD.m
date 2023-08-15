%
function [y] = oppp(n,x,p,thq);
global A

z = x;
i = 1;
while(i <= p)
        w = op(n,z);
        z = z - (1/thq(i))*w;
        i = i + 1;
end
y = x - z;
