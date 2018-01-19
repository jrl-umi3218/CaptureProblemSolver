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


%The non-linear equality constraint
%sum(delta_j/(sqrt(x(j+1))+sqrt(x(j))) - alpha sqrt(x(n)) - b
function [cineq,c,gradcineq,gradc] = constrMPC(delta, x, alpha, b, relax)
if nargin < 5, relax = false; end
n = length(x);
assert(length(delta) == n)
y = [0;x];
cineq=[];
c = -alpha*sqrt(x(n))-b;
for i=0:n-1
  c = c + delta(i+1)/(sqrt(y(i+2))+sqrt(y(i+1)));
end
if nargout > 3
  gradcineq = [];
  gradc = zeros(n,1);
  gradc(n) = -alpha/(2*sqrt(x(n)));
  for i=0:n-1
    if (i>0)
      gradc(i) = gradc(i) - delta(i+1)/(2*sqrt(y(i+1))*(sqrt(y(i+2))+sqrt(y(i+1)))^2);
    end
    gradc(i+1) = gradc(i+1) - delta(i+1)/(2*sqrt(y(i+2))*(sqrt(y(i+2))+sqrt(y(i+1)))^2);
  end
  if relax
    gradcineq = gradc;
    gradc = [];
  end
end
if relax
  cineq = c;
  c = [];
end
end