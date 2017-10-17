%remove elements with absolute value lower than threshold
function Y = clean(X,threshold)
if nargin<2
  threshold = 1e-14;
end
I = abs(X)>threshold;
Y = zeros(size(X));
Y(I) = X(I);
end