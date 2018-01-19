% Copyright 2018 CNRS-AIST JRL, CNRS-UM LIRMM
%
% This file is part of CPS.
%
% CPS is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% CPS is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Lesser General Public License for more details.
%
% You should have received a copy of the GNU Lesser General Public License
% along with CPS.  If not, see <http://www.gnu.org/licenses/>.


%Draw the problem geometry in the case n=3, in the plane phi_1 = delta(1)*g/zf;
function display3dCase(zf,zi,dzi,wimin,wimax)
g = 9.80665;
lambdamin = g/10;
lambdamax = 2*g;
s = [0, 1/3, 2/3, 1];
delta = (s(2:end).^2 - s(1:end-1).^2);
d = 1./delta;
db = d(1:end-1);
de = d(2:end);
z = zeros(length(db),1);
J = [-diag(db+de)+diag(de(1:end-1),-1) z] + [z diag(de)];
x1 = delta(1)*g/zf;
alpha = zi/g;
beta = dzi/g;

%draw in the 2d plane phi1 = x1
figure();
%  draw the boundedness constraint
phi3 = [];
r2 = 0:0.025:10;
for phi2=r2
  z = solveBoundedness(delta,alpha,beta,[x1,phi2]);
  phi3(end+1) = z;
end
plot(r2,phi3);
hold on
%  display min value for phi3
% lim = ((delta(1)/sqrt(x1)-beta)/alpha)^2;
% plot([r2(1),r2(end)],[lim,lim]);

%  draw zonotope
l = [x1;lambdamin*delta(2:end)'];
u = [x1;lambdamax*delta(2:end)'];
D = diag(u-l);
L = [1 0 0; 1 1 0; 1 1 1];
Z = L*(repmat(l,[1,5]) + D*[0 0 0 0 0; 0 0 1 1 0; 0 1 1 0 0]);
plot(Z(2,:),Z(3,:));

%  draw limits on phi3
plot([r2(1),r2(end)], [wimin^2 wimin^2])
plot([r2(1),r2(end)], [wimax^2 wimax^2])

%  draw minimum value of objective
p = -J(:,2:end)\(x1*J(:,1));
scatter(p(1),p(2),'+');

%  draw solution
x = balanceMPC(g,delta,zi,dzi,zf,lambdamin,lambdamax,wimin,wimax);
scatter(x(2),x(3),'x');

%  draw relaxed solution
x = balanceMPC(g,delta,zi,dzi,zf,lambdamin,lambdamax,wimin,wimax, true);
scatter(x(2),x(3),'o');
end