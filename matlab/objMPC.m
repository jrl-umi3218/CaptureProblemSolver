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


%Objective function of the original capture problem
%sum(((x[i+1]-x[i])\delta[i] - (x[i]-x[i-1])\delta[i-1])^2) (with x[0] = 0)
function [f,gradf] = objMPC(delta, x)
n = length(x);
assert(length(delta) == n)
y = [0;x];
f = 0;
for i=1:n-1
  f = f + ((y(i+2)-y(i+1))/delta(i+1) - (y(i+1)-y(i+0))/delta(i+0))^2;
end
if nargout > 1
  gradf = zeros(n,1);
  for j=1:n-1
    if j>1
      gradf(j-1) = gradf(j-1) - 2*(y(j+1)-y(j+0))/delta(j)^2 + 2*(y(j+2)-y(j+1))/(delta(j+1)*delta(j));
    end
    gradf(j) = gradf(j) - 2*(y(j+2)-y(j+1))/delta(j+1)^2 - 2*(y(j+2)-2*y(j+1)+y(j+0))/(delta(j+1)*delta(j)) + 2*(y(j+1)-y(j))/delta(j)^2;
    gradf(j+1) = gradf(j+1) + 2*(y(j+2)-y(j+1))/delta(j+1)^2 - 2*(y(j+1)-y(j))/(delta(j+1)*delta(j));
  end
end
end