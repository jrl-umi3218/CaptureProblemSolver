%multiply efficiently the matrix A by the nullspace basis N that 
%dedicatedNullspace(az, a1, an) would have returned.
%idx is such that N*B = B(idx,:)
function [AN,idx] = applyNullspace(A,az,a1,an)
assert(~(a1 && az(1)));
n = length(az);
na = sum(az) + a1 +an;
AN = zeros(size(A,1),n-na);
idx = zeros(n,1);
az(1) = az(1) || a1;

c = 0;
k = 1;
%skip first group of activated constraints
while k<=n && az(k)
  k = k+1;
end
for i=k:n
  if az(i)
    AN(:,c) = AN(:,c) + A(:,i);
    idx(i) = c;
  else
    c = c+1;
    if (c>n-na)
      break;
    end
    AN(:,c) = A(:,i);
    idx(i) = c;
  end
end
end