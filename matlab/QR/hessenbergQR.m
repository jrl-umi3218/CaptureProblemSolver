%perform a QR decomposition on an upper Hessenberg matrix
function [g,R] = hessenbergQR(H)
n = size(H,1);
g(n-1,1) = struct('i',0,'n',0,'c',0.,'s',0.);
R = H;
for i=1:n-1
  disp(i)
  g(i) = givensCons(R,i,i);
  R = applyGivensOnLeft(R,g(i));
end
end