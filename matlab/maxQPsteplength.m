% Given the current point x and a step p, find a maximal such that x+ap
% verify the constraints
% l_i <= x_i-x_{i-1} <= u_i for i=1..n  (x_0 = 0)
% xln <= x_n <= xun
%
% i is the number of a constraint that prevents a to be bigger.
% It is 0 if p=1, and n+1 if the blocking constraints is the last one.
% type=1 (resp. -1) if the upper (resp. lower) bound is blocking and 0 if 
% p=1. 
function [a,i,type] = maxQPsteplength(x,p,l,u,xln,xun,activeZ,activeN)
act = [activeZ;activeN];

Ax = [x(1);x(2:end)-x(1:end-1);x(end)];
Ap = [p(1);p(2:end)-p(1:end-1);p(end)];

L = ([l;xln]-Ax)./Ap;
U = ([u;xun]-Ax)./Ap;
sup = Ap>0 & ~act;
inf = Ap<0 & ~act;

alpha(sup) = U(sup);
alpha(inf) = L(inf);
alpha(act) = 1e10;
types(sup) = 1;
types(inf) = -1;

[a,i] = min(alpha);
type = types(i);
a = min(a,1);
if a==1
  i=0;
  type = 0;
end
end