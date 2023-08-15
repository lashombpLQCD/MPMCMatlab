function [corrterm] = tracecorrection(M,e_i,p,rv)

%b = sparse(size(M,1),1);
%b(1) = 1;

%[~,~,~,~,~,~,~,th] = gmresdrEIGritz(M,b,p,1,1e-2,1);
%[rv] = ModLejaComplex(th);

corrterm = 0;
for i = 1:size(e_i,2)
%[~,~,~,~,~,~,~,th] = gmresdrEIGritz(M,e_i(:,i),p,1,1e-2,1);
%[rv] = ModLejaComplex(th);
[phi] = newpoly(M,e_i(:,i),p,rv);
corrterm = corrterm + e_i(:,i)'*phi;
end

