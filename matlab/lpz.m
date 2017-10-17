%solve the lp
% min. c'x
% s.t. l_i <= x_i-x_{i-1} <= u_i for i=1..n
%      xl_j <= x_j <= xu_j for j in idx
%(by convention x_0 = 0, n = length(c))
function x = lpz(c,l,u,idx,xl,xu)
n = length(c);
assert(length(l)==n);
[A,b,lb,ub] = buildLPZ(l,u,idx,xl,xu);
x = linprog(c,A,b,[],[],lb,ub);
end