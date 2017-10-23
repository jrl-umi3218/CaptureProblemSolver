% Build a matrix with main body of the form
%  S x
%  x x x
%    x x x
%      ...
%       x x x
%         x E
% where S and E are specified by startType and endType, and the x stands
% for diagonles with terms of e
% startType:
% - 1: S=[e(1); -e(1)-e(2)]
% - 2: S=-e(1)-e(2);
% - 3: S=-e(1);
% stopType:
% - 1: E=[same;e(k)];
% - 2: E=same;
% - 3: E=[same, e(k)];
% - 4: E=-e(k)
function Jj = universalJj(e,startType,endType)
n = length(e);
e = [e(:);0];
if startType==1, s=1; else s=0; end
if endType==1, eb=1; else eb=0; end
if endType==3, er=1; else er=0; end
if startType==3, ks=1; else ks=2; end
if endType==4, ke=n; else ke=n-1; end
  
if n==0 %n==1 && endType~=4
  Jj=zeros(0,1);
  return;
else
  Jj = [zeros(s,ke-ks+2+er); 
        (diag(e(ks:ke),-1) - diag([0;e(ks:ke)+e(ks+1:ke+1)]) + diag(e(ks:ke),1)) zeros(ke-ks+2,er);
        zeros(eb,ke-ks+2+er)];
end
switch startType
  case 1, Jj(1:2,1) = [e(1);-e(1)-e(2)];
  case 2, Jj(1,1) = -e(1)-e(2);
  case 3, Jj(1,1) = -e(1);
end
if endType == 1 || endType == 3
  Jj(end,end) = e(n);
end
end