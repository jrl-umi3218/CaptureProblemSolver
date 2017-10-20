%multiply efficiently the matrix A by the nullspace basis N that 
%dedicatedNullspace(az, a1, an) would have returned.
%idx is such that N*B = B(idx,:)
function [AN,idx,activeIndex] = applyNullspace(A,az,a1,an)
assert(~(a1 && az(1)));
n = length(az);
na = sum(az) + a1 +an;
AN = zeros(size(A,1),n-na);
idx = zeros(n,1);
activeIndex = zeros(na,1);
az(1) = az(1) || a1;

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