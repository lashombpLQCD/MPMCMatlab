%
function [y] = minv(n,x);
%global D
% Diagonal preconditioning:
%y = D.\x;
%Preconditioner M
global M
% y = M\x;
% % Next line is for no preconditioning.
y = x;