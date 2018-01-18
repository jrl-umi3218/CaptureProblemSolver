%given the equation 
%sum (delta(i)/(sqrt(x(j)) + sqrt(x(j+1))) - alpha * sqrt(x(end)) = beta
%find the positive solution in x(end), given the values x(1:end-1)
function z = solveBoundedness(delta,alpha,beta,x)
n = length(delta);
a = delta(1)/sqrt(x(1)) - beta;
for i=1:n-2
  a = a + delta(i+1)/(sqrt(x(i)) + sqrt(x(i+1)));
end
b = sqrt(x(n-1));
discr = (a-alpha*b)^2 + 4*alpha*(a*b+delta(n));
if discr<0
  z = NaN;
else
  z = ((-(a-alpha*b) - sqrt(discr))/(-2*alpha))^2;
end
end