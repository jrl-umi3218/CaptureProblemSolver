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


%build the data for the feasibility problem 
% find x
% s.t. l_i <= x_i-x_{i-1} <= u_i for i=1..n
%      xl_j <= x_j <= xu_j for j in idx
function [A,b,lb,ub] = buildLPZ(l,u,idx,xl,xu)
if nargin<3
  idx=[];
  xl=[];
  xu=[];
end
n = length(l);
assert(length(u)==n);
assert(length(idx)<=n);
assert(all(idx)>0 && all(idx)<=n);
assert(length(idx)==length(xl));
assert(length(idx)==length(xu));

lb=-1e4*ones(n,1);
lb(idx) = xl;
ub=1e4*ones(n,1);
ub(idx) = xu;
A = zeros(2*n,n);
b = zeros(2*n,1);
A(1,1) = -1; b(1) = -l(1); %l(1) <= x_1
A(2,1) = 1; b(2) = u(1); % x_1 <= u(1)
for i=2:n
  % l_i <= x_i-x_{i-1}
  A(2*i-1,i-1) = 1;
  A(2*i-1,i) = -1;
  b(2*i-1) = -l(i);
  % x_i-x_{i-1} <= u_i
  A(2*i,i-1) = -1;
  A(2*i,i) = 1;
  b(2*i) = u(i);
end
end