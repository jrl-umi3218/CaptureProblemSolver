%sum(delta_j/(sqrt(x(j+1))+sqrt(x(j))) - alpha sqrt(x(n)) - b
function [cineq,c,gradcineq,gradc] = constrMPC(delta, x, alpha, b, relax)
n = length(x);
assert(length(delta) == n)
y = [0;x];
cineq=[];
c = -alpha*sqrt(x(n))-b;
for i=0:n-1
  c = c + delta(i+1)/(sqrt(y(i+2))+sqrt(y(i+1)));
end
if nargout > 3
  gradcineq = [];
  gradc = zeros(n,1);
  gradc(n) = -alpha/(2*sqrt(x(n)));
  for i=0:n-1
    if (i>0)
      gradc(i) = gradc(i) - delta(i+1)/(2*sqrt(y(i+1))*(sqrt(y(i+2))+sqrt(y(i+1)))^2);
    end
    gradc(i+1) = gradc(i+1) - delta(i+1)/(2*sqrt(y(i+2))*(sqrt(y(i+2))+sqrt(y(i+1)))^2);
  end
  if relax
    gradcineq = gradc;
    gradc = [];
  end
end
if relax
  cineq = c;
  c = [];
end
end