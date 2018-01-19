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


% return the nxn Givens rotations G such that 
% G(i:i+1,i:i+1)^T [-sqrt(i)/i;1] = [X;0]
% G is encoded as g (see g2mat)
% if i is a vector, return a vector of structures.
function g = specialGivensCons(i,n)
%we go downward only to be memory efficient: the size of g is fixed after
%the first iteration
for k=length(i):-1:1
  j = i(k);
  s.i = j;
  s.n = n;
  s.c = 1/sqrt(j+1);
  s.s = sqrt(j)*s.c;
  g(k) = s;
end
end