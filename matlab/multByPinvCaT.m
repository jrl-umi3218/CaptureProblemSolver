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


%Given the following constraints
% l_i <= x_i-x_{i-1} <= u_i for i=1..n  (x_0 = 0)
% xln <= x_n <= xun
%we note Ca the matrix corresponding to the activated constraints,
%specified by  az and an.
%The function computes Y = pinv(Ca^T)*X in an optimized way.
function Y = multByPinvCaT(X,az,an)
n = length(az);
assert(size(X,1) == n);
na = sum(az) + an;
assert(na<=n); %we do not treat the case where all constraints are active.
Y = zeros(na,size(X,2));

if az(1)
  k=1;
  while k<na && az(k+1)
    k = k+1;
  end
  Y(1:k,:) = applyFirst(X,k);
  k0 = k+1;
else
  k0 = 1;
end
if an
  k=1;
  while k<na && az(n+1-k)
    k=k+1;
  end
  Y(na-k+1:na,:) = applyLast(X,k);
  kn = n-k;
else
  kn = n;
end
k=0;
for i=k0:kn
  if az(i)
    k = k+1;
  else
    if k>0
      Y(k0:k0+k-1,:) = apply(X,i-k-2,k+1);
      k0 = k0+k;
      k=0;
    end
  end
end
if k>0
  Y(k0:k0+k-1,:) = apply(X,kn-k-1,k+1);
end
end

function Y = applyFirst(X,k)
Y = zeros(k,size(X,2));
Y(k,:) = X(k,:);
for i=k-1:-1:1
  Y(i,:) = Y(i+1,:) + X(i,:);
end
end

function Y = applyLast(X,k)
Y = zeros(k,size(X,2));
if k==1
  Y = X(end,:);
else
  n = size(X,1);
  Y(1,:) = -X(n-k+1,:);
  for i=2:k-1
    Y(i,:) = Y(i-1,:) - X(n-k+i,:);
  end
  Y(k,:) = -Y(k-1,:) + X(end,:);
end

end

%apply to X(s+1:k,:)
function Y = apply(X,s,k)
Y = zeros(k-1,size(X,2));
for i=1:k-1
  Y(i,:) = [-(k-i)*ones(1,i) i*ones(1,k-i)]*X(s+(1:k),:);
end
Y = Y/k;
end