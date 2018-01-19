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


%find a point x in the polytope
% l_i <= x_i-x_{i-1} <= u_i for i=1..n
% xln <= x_n <= xun
%
%b indicates if the constraints are feasible
%
%The points verifying the first n equalities are exactly
%x = L*l + sum(lambda_i*(u_i-l_i)*L_i) with lambda_i in [0,1], L the lower
%triangular matrix filled by ones, and L_i its i-th column.
%This describes a zonotope Z.
%
%We thus have x_n = sum(l_i) + sum(lambda_i*(u_i-l_i)).
%Let x(a) = L*l + sum(a*(u_i-l_i)*L_i). For a in [0,1], x(a) is in Z
%We have x_n(a) =  sum(l_i) + sum(a*(u_i-l_i)) = sum(l_i) + a*sum(u_i-l_i)
%x(0) is the point of Z with the lowest x_n, while x(1) has the highest 
%x_n.
%We compute al and au such that x_n(al) = xln and x_n(au) = xun
%The intersection of Z with the last constraint is non empty if [0,1]
%intersects [al,au]. Let a0 be the value in the middle of this
%intersection, we return x(a0).

function [x,b] = initialPoint(l,u,xln,xun)
assert(xln<=xun);
assert(all(l<=u));
s = sum(l);
d = u(:)-l(:);
t = sum(d);
if abs(t) < 1e-15
  b = xln<=s && s<=xun;
  x = l(:);
else
  al = (xln-s)/t;
  au = (xun-s)/t;
  b = al<=1 && au>=0;
  L = tril(ones(length(l)));
  x = L*(l(:) + (max(al,0)+min(au,1))/2*d);
end
end