%objective of the balance problem
function [f,gradf] = objMPC(delta, x)
n = length(x);
assert(length(delta) == n)
y = [0;x];
f = 0;
for i=1:n-1
  f = f + ((y(i+2)-y(i+1))/delta(i+1) - (y(i+1)-y(i+0))/delta(i+0))^2;
end
if nargout > 1
  gradf = zeros(n,1);
  for j=1:n-1
    if j>1
      gradf(j-1) = gradf(j-1) - 2*(y(j+1)-y(j+0))/delta(j)^2 + 2*(y(j+2)-y(j+1))/(delta(j+1)*delta(j));
    end
    gradf(j) = gradf(j) - 2*(y(j+2)-y(j+1))/delta(j+1)^2 - 2*(y(j+2)-2*y(j+1)+y(j+0))/(delta(j+1)*delta(j)) + 2*(y(j+1)-y(j))/delta(j)^2;
    gradf(j+1) = gradf(j+1) + 2*(y(j+2)-y(j+1))/delta(j+1)^2 - 2*(y(j+1)-y(j))/(delta(j+1)*delta(j));
  end
end
end