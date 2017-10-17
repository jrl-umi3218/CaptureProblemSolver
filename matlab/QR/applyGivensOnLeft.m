% perform Y = G_k^T G_{k-1}^T... G_1^T X, where G_i is the Givens rotation 
% described by g(i) (see g2mat).
function Y = applyGivensOnLeft(X,g)
Y = X;
for i=1:length(g)
  gi = g(i);
  Y(gi.i+(0:1),:) = [gi.c -gi.s; gi.s gi.c]* Y(gi.i+(0:1),:);
end
end