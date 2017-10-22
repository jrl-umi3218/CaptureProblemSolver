function JN = directJN(delta,az,a1,an)
assert(~(a1 && az(1)));
n = length(az);
na = sum(az) + a1 +an;
JN = zeros(n-1,n-na);
az(1) = az(1) || a1;
az(end+1) = an;
d = 1./delta;

if na==n, return; end
if na==0, JN = universalJj(d,2,3); return; end

a0=0;
while az(a0+1) && a0<=n
  a0 = a0+1;
end
k = a0;
i1=0;
while k+i1<=n && ~az(k+i1+1)
  i1 = i1+1;
end
if a0>0
  JN(k:k+i1-2,1:i1-1) = buildJ0(d(k+1:k+i1-1), 0);
elseif i1==n
  JN = universalJj(d, 2, 2);
elseif i1==(n+1-na)
  %JN = [buildJp(d(1:n+1-na), 2, true); zeros(na-2,n-na)];
  JN(1:n+1-na,:) = universalJj(d(1:n+1-na),2,1);
else
  JN(1:i1,1:i1) = universalJj(d(1:i1),2,4);
end
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


function J0 = buildJ0(e, c)
if c==0
  if length(e)==1
    J0 = [e;-e];
  else
    J0 = [diag(e) - diag(e(1:end-1)+e(2:end),-1) + diag(e(2:end-1),-2); 
          zeros(1, length(e)-2) e(end) -e(end)];
  end
elseif c==1
  J0 = -diag(e(1:end-1)+e(2:end)) + diag(e(2:end-1),-1) + diag(e(2:end-1),1);
else
  J0 = -diag(e(1:end) + [e(2:end); 0]) + diag(e(2:end),1) + diag(e(2:end),-1);
end
end

function Jj = buildJj(e)
if length(e)==1
  Jj = [-e e; e -e];
else
  Jj = -diag([0;e]+[e;0]) + diag(e,-1) + diag(e, 1);
end
end

function Jp = buildJp(e, c, s)
if s==true
  e0 = e(1);
  e = e(2:end);
end
n = length(e);
if c == 0
  Jp = zeros(n,n+1);
  if n>=1
    Jp(1,1:2) = [-e(1),e(1)];
    Jp(2:end,:) = [diag(e(1:end-1)) zeros(n-1,2)]... 
                - [zeros(n-1,1) diag(e(1:end-1)+e(2:end)) zeros(n-1,1)]... 
                + [zeros(n-1,2) diag(e(2:end))];
  end
elseif c == 1
  Jp = -diag([0;e(1:end-1)] + e(1:end)) + diag(e(1:end-1),-1) + diag(e(1:end-1),1);
else
  Jp = zeros(n+1,n);
  Jp(1:2,1) = [-e(1); e(1)];
  Jp(:,2:end) = [diag(e(1:end-1)); zeros(2,n-1)]...
              - [zeros(1,n-1); diag(e(1:end-1) + e(2:end)); zeros(1,n-1)]...
              + [zeros(2,n-1); diag(e(2:end))];
end
if s==true
  Jp(1,1) = Jp(1,1)-e0;
end
end