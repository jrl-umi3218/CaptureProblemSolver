%solve the problem 
% min. ||c+j'x||^2
% s.t. l_i <= x_i-x_{i-1} <= u_i for i=1..n  (x_0 = 0)
%      xln <= x_n <= xun
%
function [x,lambda] = feasibilityLS(j,c,l,u,xln,xun)
n = length(j);
assert(length(l) == n);
assert(length(u) == n);
assert(all(l<=1e-14) && all(-1e-14<=u)); %x=0 is supposed to be a feasible point
assert(xln<=0 && 0<=xun);

%zonotope constraints: -1 when active at lower bound, +1 if active at upper
%bound, 2 if active as equality, 0 if inactive;
activeZ = zeros(n,1);
activeN = 0; %constraint on x_n: -1 if active at lower bound, +1...
activeZ(abs(l-u)<=1e-14)=2;
activeN(xln==xun)=2;

lambda = zeros(n+1,1);

x = zeros(n,1);
%finished = checkKKT(j,c,x,lambdaZ,lambdaN,l,u,xln,xun,1e-6,1e-6);
maxIter = 10*n;
for k=1:maxIter
  [jn,idx] = applyNullspace(j',activeZ~=0,false,activeN~=0);
  %z = -pinv(jn)*(c+j'*x);
  z = -(jn'/(jn*jn'))*(c+j'*x);
  %compute p = N*z
  p = zeros(n,1);
  nz = idx~=0;
  p(nz) = z(idx(nz));
  if norm(p,Inf)<1e-12
    act = [activeZ;activeN];
    %[~,Cact] = dedicatedNullspace(activeZ~=0,false,activeN~=0);
    %lambda(act~=0) = -(j'*x+c)*(pinv(Cact')*j);
    lambda(act~=0) = -(j'*x+c)*(multByPinvCaT(j,activeZ~=0,false,activeN~=0));
    lambda(act==0) = 0;
    if all(lambda(act==-1)<=0) && all(lambda(act==1)>=0)
      return;
    else
      %deactivate
      tmp = zeros(n+1,1);
      tmp(act==-1) = -lambda(act==-1);
      tmp(act== 1) = lambda(act==1);
      [~,i] = min(tmp);
      if (i<=n)
        activeZ(i) = 0;
      else
        activeN = 0;
      end
    end
  else
    [a,i,t] = maxQPsteplength(x,p,l,u,xln,xun,activeZ~=0,activeN~=0);
    x = x + a*p;
    if i>0
      if i<=n
        activeZ(i) = t;
      else
        activeN = t;
      end
    end
  end
end
end

