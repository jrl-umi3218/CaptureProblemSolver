function x = unconstrainedGN(x0,delta,alpha,b,w)
n = length(x0);
delta = delta(:);
%Form the matrix A such that the objective is ||Ax||^2
A = [diag(-1./delta(2:end)-1./delta(1:end-1)) zeros(n-1,1)]... 
  + [zeros(n-1,1) diag(1./delta(2:end))]...
  + [diag(1./delta(2:end-1),-1) zeros(n-1,1)];

%Form the QR factorization of A
[Qa,Ra] = qr(A); %TODO optimize, knowing that A is tridiagonal

[f,g] = constrMPC(delta,x0,alpha,b);

end