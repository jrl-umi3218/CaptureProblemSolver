function b = evalBoundedness(x,delta,alpha,beta)
n = length(x);
assert(length(delta) == n);
b = delta(1)/sqrt(x(1)) - beta;
for i=1:n-1
  b = b + delta(i+1)/(sqrt(x(i))+sqrt(x(i+1)));
end
b = b - alpha*sqrt(x(end));
end