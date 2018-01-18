%Build and solve the capture problem
%If relax is true (false by default), build a relaxed version in which the
%non-linear equlity constraint c=0 is turned into the convex constraint
%c<=0%
%
%Example
% delta = [0.01, 0.03, 0.05, 0.07, 0.09, 0.11, 0.13, 0.15, 0.17, 0.19];
% g = 9.80665;
% lambda_max = 19.6133;
% lambda_min = 0.980665;
% omega_i_max = 3.5004454562949268;
% omega_i_min = 2.7492602295351802;
% z_bar = 0.91307774471737957;
% z_f = 0.8;
% zd_bar = -0.17413705390465162;
% balanceMPC(g,delta,z_bar,zd_bar,z_f,lambda_min,lambda_max, omega_i_min, omega_i_max)
function phi = balanceMPC(g,delta,zi,dzi,zf,lmin,lmax,wimin,wimax,relax)
assert(nargin>=8)
twoD = nargin==8;
delta=delta(:);
n = length(delta);

if twoD
  alpha = 0;
  b = (dzi+wimin*zi)/g;
else
  alpha = zi/g;
  b = dzi/g;
end

if nargin<10
  relax = false;
end

o = @(x) objMPC(delta,x);
c = @(x) constrMPC(delta,x,alpha,b,relax);

A = diag(ones(n,1))+diag(-ones(n-1,1),-1);
A = [-A;A];
b = [-lmin*delta; lmax*delta];
xl = zeros(n,1);
xu = 1000*ones(n,1);
if twoD
  Aeq = zeros(2,n);
  Aeq(1,1) = 1; Aeq(2,n) = 1;
  beq = [delta(1)*g/zf;wimin^2];
else
  Aeq = zeros(1,n);
  Aeq(1,1) = 1;
  beq = delta(1)*g/zf;
  xl(n) = wimin^2;
  xu(n) = wimax^2;
end
L = tril(ones(n,n));
phi0 = (lmin+lmax)*L*delta/2;
opt = optimset('GradObj','on','GradConstr','on');
phi = fmincon(o,phi0,A,b,Aeq,beq,xl,xu,c,opt);

end
