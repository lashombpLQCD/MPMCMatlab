
% z4noise2 takes a random vector on (0,1) and converts it to Z4 noise.

% z2noise.m
% Generates a Z2 noise vector of size n
% Assumptions: n is given as a positive integer
%
% Input:  {n,m}>1    length of the Z2 noise vector
% Output: b      Z2 noise vector

function [b] = z4noise(n,b)

% rand('seed',0)
% b=rand(n,m);

  for i=1:n
    if (b(i,1) > .75)
      b(i,1) =  1;
    elseif(b(i,1) > .5)
      b(i,1) = -1;
    elseif(b(i,1) > .25)
      b(i,1) = 1i;
    else
      b(i,1) = -1i;
    end
  end
