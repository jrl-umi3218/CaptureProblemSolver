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