%Considering that G is defined with respect to the submatrix S of a matrix 
%A, return G' such that G'^T A has the same effect on S than G^T S.
%A has n rows, and S begin at row i
function g = expandGivens(g,i,n)
for k=1:length(g)
  g(k).i = g(k).i+i-1;
  g(k).n = n;
end
end