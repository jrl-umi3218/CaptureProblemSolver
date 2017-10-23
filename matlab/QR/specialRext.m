function R = specialRext(e,endType)
n = length(e);
d = -(1:n)';
l = 1./sqrt((1:n).*(2:n+1))';
if endType==2 || endType==3
  d(end) = 0;
  l(end) = 1/sqrt(n);
end
c1 = l.*(d-1).*e;
c2 = l.*d.*[e(2:end);0];
switch endType
  case 1
    R = [diag(c1) - diag(c1(1:end-1)+c2(1:end-1),1) + diag(c2(1:end-2),2);zeros(1,n)];
  case 2
    R = diag(c1) - diag(c1(1:end-1)+c2(1:end-1),1) + diag(c2(1:end-2),2);
  case 3
    R = [diag(c1) zeros(n,1)] - [zeros(n,1) diag(c1+c2)] + [zeros(n,1) diag(c2(1:end-1),1)];
  case 4
    R = [[diag(c1) zeros(n,1)] - [zeros(n,1) diag(c1+c2)] + [zeros(n,1) diag(c2(1:end-1),1)]; zeros(1,n+1)]; 
end
end