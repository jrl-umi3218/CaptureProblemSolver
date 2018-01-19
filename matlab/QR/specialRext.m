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


%Extension of specialR to more types of ending
%See universalJj for the meaning of endType
function R = specialRext(e,endType)
n = length(e);
d = -(1:n)';
l = 1./sqrt((1:n).*(2:n+1))';
if endType==2 || endType==3
  d(end) = 0;
  l(end) = 1/sqrt(n);
end
c1 = l.*(d-1).*e;
c2 = l.*d.*[e(2:end);0];
switch endType
  case 1
    R = [diag(c1) - diag(c1(1:end-1)+c2(1:end-1),1) + diag(c2(1:end-2),2);zeros(1,n)];
  case 2
    R = diag(c1) - diag(c1(1:end-1)+c2(1:end-1),1) + diag(c2(1:end-2),2);
  case 3
    R = [diag(c1) zeros(n,1)] - [zeros(n,1) diag(c1+c2)] + [zeros(n,1) diag(c2(1:end-1),1)];
  case 4
    R = [[diag(c1) zeros(n,1)] - [zeros(n,1) diag(c1+c2)] + [zeros(n,1) diag(c2(1:end-1),1)]; zeros(1,n+1)]; 
end
end