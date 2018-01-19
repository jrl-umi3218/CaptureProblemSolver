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


%solve the lp
% min. c'x
% s.t. l_i <= x_i-x_{i-1} <= u_i for i=1..n
%      xl_j <= x_j <= xu_j for j in idx
%(by convention x_0 = 0, n = length(c))
function x = lpz(c,l,u,idx,xl,xu)
n = length(c);
assert(length(l)==n);
[A,b,lb,ub] = buildLPZ(l,u,idx,xl,xu);
x = linprog(c,A,b,[],[],lb,ub);
end