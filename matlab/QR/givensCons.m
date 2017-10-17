%return the Givens rotation zeroing A(i+1,j) with A(i,j)
%g is a structure (i,n,c,s)
%g described a nxn matrix that can be obtained with
%G = eye(n); G(i:i+1,i:i+1) = [c s; -s c]
function g = givensCons(A,i,j)
g.i = i;
g.n = size(A,1);
a = A(i,j);
b = A(i+1,j);
if b==0
  g.c = 1;
  g.s = 0;
else
  if abs(b)>abs(a)
    t = -a/b;
    s = 1/sqrt(1+t^2);
    g.c = s*t;
    g.s = s;
  else
    t = -b/a;
    g.c = 1/sqrt(1+t^2);
    g.s = g.c*t;
  end
end
end