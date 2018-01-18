%solve the problem
% min. ||c(x)||^2
% s.t. l_i <= x_i-x_{i-1} <= u_i for i=1..n  (x_0 = 0)
%      xln <= x_n <= xun
%
%where c is the boundedness constraint
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
% feasibilitySQP(g,delta,z_bar,zd_bar,z_f,lambda_min,lambda_max, omega_i_min, omega_i_max)
function x = feasibilitySQP(g,delta,zi,dzi,zf,lmin,lmax,wimin,wimax)
n = length(delta);
l = lmin*delta(:);
u = lmax*delta(:);
l(1) = delta(1)*g/zf;
u(1) = l(1);
xln = wimin^2;
xun = wimax^2;
alpha = zi/g;
b = dzi/g;
c = @(x) constrMPC(delta,x,alpha,b);
[x,b] = initialPoint(l,u,xln,xun);
assert(b, 'could not find an initial point');
lambda = zeros(n+1,1);

maxIter = 100;
for i=1:maxIter
  [~,f,~,g] = c(x);
  disp(['||c||^2 : ' num2str(f^2)])
  if checkKKT(f,g,x,lambda,l,u,xln,xun,1e-6,1e-6)
%    disp(norm(f));
    break;
  end
  %solve the LS
  Ax = [x(1); x(2:end)-x(1:end-1)];
  [p,pl] = feasibilityLS(g,f,l-Ax,u-Ax,xln-x(n),xun-x(n));
  
  %line search
  a = 1;
  beta = 0.9;
  c1 = 0.01;
  cgp = c1*f*g'*p;
  while 0.5*c(x+a*p)^2 > 0.5*f^2 + a*cgp
    a = beta*a;
    if a<1e-8
      error('step is too small');
    end
  end
  x = x+a*p;
  lambda = lambda + a*(lambda-pl);
end
disp(norm(f));
end


function b = checkKKT(f,g,x,lambda,l,u,xln,xun,tp,td)
%see Stan's thesis 4.3.5
tx = tp*(1+norm(x,Inf));
tl = td*(1+norm(lambda,Inf));

en = zeros(length(l),1); en(end) = 1; %n-th column of the nxn identity
gradL = f*g + [lambda(1:end-2)-lambda(2:end-1);lambda(end-1)] + lambda(end)*en;

cstr = [x(1); x(2:end)-x(1:end-1);x(end)];
kktcstr = all(   (abs(cstr-[l;xln])<=tx & lambda<=tl) ...
               | ([l;xln]-tx<=cstr & cstr<=[u;xun]+tx & abs(lambda)<=tl) ...
               | (abs(cstr-[u;xun])<=tx & lambda>=-tl));
b = norm(gradL,Inf) <= tl && kktcstr; 
end