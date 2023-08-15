%
function [Y] = opppp(n,pp,thq);
global A

% Find p(A) for matrix A.

% New version avoids using coeff's.   Very slick.

% Newton:  (actually uses roots 0, theta_1, theta_2, ... theta_deg-2, where 
%             the degree is deg-1 since it is poly p, not q, don't see 0 )

% y = cp(1)*x;
% pr = x;
% for j=2:pp
%     pr =  op(n,pr)-thq(j-1)*pr;
%     y = y + cp(j)*pr;
% end

I = speye(n,n);
Y = (1/thq(1))*I;
W = Y;
for j=2:pp
    W = (1/thq(j))*(thq(j-1)*W - A*W);
    Y = Y + W;                              
end





