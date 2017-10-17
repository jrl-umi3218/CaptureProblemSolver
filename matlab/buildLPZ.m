%build the data for the feasibility problem 
% find x
% s.t. l_i <= x_i-x_{i-1} <= u_i for i=1..n
%      xl_j <= x_j <= xu_j for j in idx
function [A,b,lb,ub] = buildLPZ(l,u,idx,xl,xu)
if nargin<3
  idx=[];
  xl=[];
  xu=[];
end
n = length(l);
assert(length(u)==n);
assert(length(idx)<=n);
assert(all(idx)>0 && all(idx)<=n);
assert(length(idx)==length(xl));
assert(length(idx)==length(xu));

lb=-1e4*ones(n,1);
lb(idx) = xl;
ub=1e4*ones(n,1);
ub(idx) = xu;
A = zeros(2*n,n);
b = zeros(2*n,1);
A(1,1) = -1; b(1) = -l(1); %l(1) <= x_1
A(2,1) = 1; b(2) = u(1); % x_1 <= u(1)
for i=2:n
  % l_i <= x_i-x_{i-1}
  A(2*i-1,i-1) = 1;
  A(2*i-1,i) = -1;
  b(2*i-1) = -l(i);
  % x_i-x_{i-1} <= u_i
  A(2*i,i-1) = -1;
  A(2*i,i) = 1;
  b(2*i) = u(i);
end
end