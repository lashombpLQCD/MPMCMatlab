%
function [y] = opppp(n,x,pp,p2,thq,thq2);
global A

% SADLY NOT done yet for keeping real arith in the complex case.

% New version avoids using coeff's.   Very slick.


y = (1/thq2(1))*x;
w = y;
for j=2:p2
    w = (1/thq2(j))*(thq2(j-1)*w - oppprQCD(n,w,pp,thq));
    y = y + w;                              
end





