% return the nxn Givens rotations G such that 
% G(i:i+1,i:i+1)^T [-sqrt(i)/i;1] = [X;0]
% G is encoded as g (see g2mat)
% if i is a vector, return a vector of structures.
function g = specialGivensCons(i,n)
%we go downward only to be memory efficient: the size of g is fixed after
%the first iteration
for k=length(i):-1:1
  j = i(k);
  s.i = j;
  s.n = n;
  s.c = 1/sqrt(j+1);
  s.s = sqrt(j)*s.c;
  g(k) = s;
end
end