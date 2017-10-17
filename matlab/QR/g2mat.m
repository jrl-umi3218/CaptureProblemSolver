%g = (i,n,c,s) describes a nxn matrix that can be obtained with
%G = eye(n); G(i:i+1,i:i+1) = [c s; -s c]
function G = g2mat(g)
G = eye(g.n);
r = g.i + (0:1);
G(r,r) = [g.c g.s; -g.s g.c];
end