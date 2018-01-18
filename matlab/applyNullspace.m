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