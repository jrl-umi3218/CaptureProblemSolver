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


%return the Givens rotation zeroing A(i+1,j) with A(i,j)
%g is a structure (i,n,c,s)
%g described a nxn matrix that can be obtained with
%G = eye(n); G(i:i+1,i:i+1) = [c s; -s c]
function g = givensCons(A,i,j)
g.i = i;
g.n = size(A,1);
a = A(i,j);
b = A(i+1,j);
if b==0
  g.c = 1;
  g.s = 0;
else
  if abs(b)>abs(a)
    t = -a/b;
    s = 1/sqrt(1+t^2);
    g.c = s*t;
    g.s = s;
  else
    t = -b/a;
    g.c = 1/sqrt(1+t^2);
    g.s = g.c*t;
  end
end
end