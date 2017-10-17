%solve the problem 
% min. ||Ax-b||^2 + w^2*||c+j'x||^2
% s.t. l_i <= x_i-x_{i-1} <= u_i for i=1..n  (x_0 = 0)
%      xln <= x_n <= xun
%
% A is supposed to be (n-1) by n with the particular structure:
%    | -a0-a1    a1    0     0  ...   0       0        0 |
%    |   a1   -a1-a2   a2    0  ...   0       0        0 |
%    |   0       a2  -a2-a3  a3       0       0        0 |
%    |                ...                                |
%    |                             a{n-1} -a{n-1}-an  an |
%
% x = 0 is supposed to be feasible for the inequalities
function x = dedicatedLS(A,b,j,c,w,l,u,xln,xun)
n = size(A,2);
assert(size(A,1) == n-1);
assert(length(b) == n-1);
assert(length(j) == n);
assert(length(l) == n);
assert(length(u) == n);

%zonotope constraints: -1 when active at lower bound, +1 if active at upper
%bound, 0 if inactive;
activeZ = zeros(n,1);
activeN = 0; %constraint on x_n: -1 if active at lower bound, +1...

finished = false;
while ~finished
  
end
end