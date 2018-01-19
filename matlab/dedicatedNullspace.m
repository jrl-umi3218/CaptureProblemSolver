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


% Compute the nullspace of active constraints among
% x_i-x_{i-1} = a_i for i=1..n  (x_0 = 0) (Z)
% x_1 = x1                                (1)
% x_n = xn                                (n)
%
% az is a n-vector of booleans specifying which constraints of Z are
% active.
% a1 is a boolean for the activity of (1), and an for the activity of (n)
%
% N is the null space basis
% Aa is the matrix of active constraints
function [N,Aa] = dedicatedNullspace(az, a1, an)
assert(~(a1 && az(1)));
n = length(az);
na = sum(az) + a1 +an;
N = zeros(n, n-na);

az(1) = az(1) || a1;
if (any(az) || an)
  p = processAZ(az);
  if p(1,1)==1
    c = p(1,2)+1;
    i0=2;
  else
    c=1;
    i0=1;
  end
  r=1;
  for i=i0:size(p,1)
    %free variables
    l = p(i,1)-2 - c;
    if l>=0
      N(c+[0:l],r+[0:l]) = eye(l+1);
      r = r + l+1;
    end
    %for j=c:p(i,1)-2
    %  N(j,r) = 1;
    %  r = r+1;
    %end
    c = p(i,1)-1;
    if c+p(i,2)==n && an
      c=n+1;
      break;
    end
    N(c:c+p(i,2),r) = 1;
    c = p(i,1)+ p(i,2);
    r = r+1;
  end
  for j=c:n-1
    N(j,r) = 1;
    r = r+1;
  end
  if ~an && c<=n
    N(n,r)=1;
  end
else
  N = eye(n);
end
  
if nargout>1
  Aa = zeros(na,n);
  r = 1;
  if a1 || az(1)
    Aa(1,1) = 1;
    r=r+1;
  end
  for i=2:n
    if az(i)
      Aa(r,i-1:i) = [-1 1];
      r = r+1;
    end
  end
  if an
    Aa(r,n) = 1;
  end
end
end

function p=processAZ(az)
p = zeros(0,2);
k=0;
past = false;
for i=1:length(az)
  if az(i)
    if past
      p(k,2) = p(k,2)+1;
    else
      k = k+1;
      past = true;
      p(k,:) = [i,1];
    end
  else
    past = false;
  end
end
end