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


%Compute the triangular matrix of the QR decomposition of A, where A has
%the form
% | -e0    e0                                      |
% |  e0  -e0-e1   e1                               |
% |        e1   -e1-e2  e2                         |
% |               ...                              |
% |                     e_{k-1}  -e_{k-1}-e_k  e_k |
% |                                  e_k      -e_k |
% when ext is false (default) and the last element -e_k is replaced by
% -e_k-e_{k+1} if ext is true
%
% The matrix R is such that Q*R = A where Q is obtained with
% specialGivensCons
function R = specialR(A,ext)
if nargin<2
  ext = false;
end
e = [diag(A,-1);0];
if ext
 e(end) = -A(end,end)-A(end-1,end);
end
d = -(1:length(e)-1)';
l = 1./sqrt((1:length(e)-1).*(2:length(e)))';
c1 = l.*(d-1).*e(1:end-1);
c2 = l.*d.*e(2:end);
R = diag([c1;0])-diag(c1+c2,1) + diag(c2(1:end-1),2);
if ext
 R(end) = (A(end,end)+A(end-1,end))/sqrt(length(e));
end
end