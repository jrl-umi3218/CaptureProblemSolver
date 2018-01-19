% Copyright 2018 CNRS-AIST JRL, CNRS-UM LIRMM
%
% This file is part of CPS.
%
% CPS is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% CPS is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Lesser General Public License for more details.
%
% You should have received a copy of the GNU Lesser General Public License
% along with CPS.  If not, see <http://www.gnu.org/licenses/>.


% Given the current point x and a step p, find a maximal such that x+ap
% verify the constraints
% l_i <= x_i-x_{i-1} <= u_i for i=1..n  (x_0 = 0)
% xln <= x_n <= xun
%
% i is the number of a constraint that prevents a to be bigger.
% It is 0 if p=1, and n+1 if the blocking constraints is the last one.
% type=1 (resp. -1) if the upper (resp. lower) bound is blocking and 0 if 
% p=1. 
function [a,i,type] = maxQPsteplength(x,p,l,u,xln,xun,activeZ,activeN)
act = [activeZ;activeN];

Ax = [x(1);x(2:end)-x(1:end-1);x(end)];
Ap = [p(1);p(2:end)-p(1:end-1);p(end)];

L = ([l;xln]-Ax)./Ap;
U = ([u;xun]-Ax)./Ap;
sup = Ap>0 & ~act;
inf = Ap<0 & ~act;

alpha(sup) = U(sup);
alpha(inf) = L(inf);
alpha(act) = 1e10;
types(sup) = 1;
types(inf) = -1;

[a,i] = min(alpha);
type = types(i);
a = min(a,1);
if a==1
  i=0;
  type = 0;
end
end