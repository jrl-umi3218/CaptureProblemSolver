%Build J*N directly (no multiplication)
%
%Example
% delta = [0.01, 0.03, 0.05, 0.07, 0.09, 0.11, 0.13, 0.15, 0.17, 0.19];
% directJN(Delta,[1 1 0 0 1 0 0 1 1 1],false)
function JN = directJN(delta,az,an)
n = length(az);
na = sum(az) + an;
JN = zeros(n-1,n-na);
az(end+1) = an;
d = 1./delta;

if na==n, return; end
%if na==0, JN = universalJj(d,2,3); return; end

%number of consecutive activated constraints at the end
ap=0;
while az(n+1-ap)
  ap = ap+1;
end

a0=0;
while az(a0+1) && a0<=n
  a0 = a0+1;
end
if a0>0
  %JN(k:k+i1-2,1:i1-1) = buildJ0(d(k+1:k+i1-1), 0);
  startType = 1;
  startOffset = a0-1;
  addDim = 1;
else
  startType = 2;
  startOffset = 0;
  addDim = 0;
end

if a0+ap==na
  if ap==0
    if na~=n-1
      JN(1+startOffset:end,:) = universalJj(d(1+a0:end), startType, 3);
    else
      JN(1+startOffset:end,:) = d(end);
    end
  elseif ap==1
    JN(1+startOffset:end,:) = universalJj(d(1+a0:end), startType, 2);
  else
    JN(1+startOffset:end+2-ap,:) = universalJj(d(1+a0:end+1-ap), startType, 1);
  end
  return;
end

k = a0;
i1=0;
while k+i1<=n && ~az(k+i1+1)
  i1 = i1+1;
end
JN(startOffset+(1:i1+addDim),1:i1) = universalJj(d(a0+(1:i1)),startType,4);

k = k + i1;
c = i1;
a1 = 0;
while k+a1<=n && az(k+a1+1) 
  a1 = a1+1;
end
k = k + a1;
while true
  i2=0;
  while k+i2<=n&& ~az(k+i2+1) 
    i2 = i2+1;
  end
  a2 = 0;
  while k+i2+a2<=n && az(k+i2+a2+1) 
    a2 = a2+1;
  end
  if (k+i2+a2==n+1)
    if a2==0
      JN(k:k+i2-2,c:c+i2-1) = universalJj(d(end-(i2-2):end),3,3);
    elseif a2==1
      JN(k:k+i2-1,c:c+i2-1) = universalJj(d(end-i2+1:end),3,2);
    else
      JN(k:k+i2,c:c+i2-1) = universalJj(d(end-i2-a2+2:end-a2+1),3,1);
    end
    break;
  end
  JN(k:k+i2,c:c+i2) = universalJj(d(k+1:k+i2),3,4);
  c = c+i2;
  k = k + i2 + a2;
end
end
