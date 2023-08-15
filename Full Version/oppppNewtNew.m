%
function [y] = opppp(n,x,pp,thq);
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

y = (1/thq(1))*x;
w = y;
for j=2:pp
%     wv = op(n,w);  %not sure what this is, maybe started to change it?
%     wv2 = 
    w = (1/thq(j))*(thq(j-1)*w - op(n,w));
    y = y + w;                              
end


    
% %Regular power basis:
% co = th;
% y = co(1)*x;
% z = x;
% for i=2:p
%     z = op(n,z);
%     y = y + co(i)*z;
% end



