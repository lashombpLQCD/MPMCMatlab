%
function [p] = opppp(n,x,deg,thq);
global A

% New version avoids using coeff's.   Very slick.

% Newton:  (actually uses roots 0, theta_1, theta_2, ... theta_deg-2, where 
%             the degree is deg-1 since it is poly p, not q, don't see 0 )

% y = cp(1)*x;
% pr = x;
% for j=2:pp
%     pr =  op(n,pr)-thq(j-1)*pr;
%     y = y + cp(j)*pr;
% end

prod = x;
p = zeros(n,1);
j = 1;

while( j <= deg-1 )
    threcip = 1/thq(j);
    ttp = threcip*prod;
    p = p + ttp;
    prod = prod - op(n,ttp);
    j = j + 1;
end
    threcip = 1/thq(deg);
    p = p + threcip*prod;


