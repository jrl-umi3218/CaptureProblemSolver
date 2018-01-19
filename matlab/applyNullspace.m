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


%multiply efficiently the matrix A by the nullspace basis N that 
%dedicatedNullspace(az, an) would have returned (see dedicatedNullspace for
%more details).
%idx is such that N*B = B(idx,:)
%
%Example
% MN = applyNullspace([1 1 0 0 1 0 0 1 1 1],false)
function [AN,idx,activeIndex] = applyNullspace(A,az,an)
n = length(az);
na = sum(az) +an;
AN = zeros(size(A,1),n-na);
idx = zeros(n,1);
activeIndex = zeros(na,1);

c = 0;
k = 1;
%skip first group of activated constraints
while k<=n && az(k)
  activeIndex(k) = k;
  k = k+1;
end
ca = k;
for i=k:n
  if az(i)
    AN(:,c) = AN(:,c) + A(:,i);
    idx(i) = c;
    activeIndex(ca) = i;
    ca = ca+1;
  else
    c = c+1;
    if (c>n-na)
      break;
    end
    AN(:,c) = A(:,i);
    idx(i) = c;
  end
end
activeIndex(ca:end) = n+ca+1-na:n+1;