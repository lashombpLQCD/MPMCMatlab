%
function [y] = oppp(n,x,p,thq);
global A

z = x;
i = 1;
while(i <= p)
    if( imag(thq(i)) == 0 )
        w = op(n,z);
        z = z - (1/thq(i))*w;
        i = i + 1;
    else
        alp = real(thq(i));
        bet = imag(thq(i));
        mod = alp*alp + bet*bet;
        w = op(n,z);
        w2 = op(n,w);
        temp = w2 - 2*alp*w;
        z = (1/mod)*temp + z;
        i = i + 2;
    end
end
y = x - z;
