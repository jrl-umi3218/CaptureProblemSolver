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


%Build the objective matrix of the original problem
%
%Example
% delta = [0.01, 0.03, 0.05, 0.07, 0.09, 0.11, 0.13, 0.15, 0.17, 0.19];
% J = buildObjMatrix(delta)
function J = buildObjMatrix(delta)
n = length(delta);
J = [diag(-1./delta(2:end)-1./delta(1:end-1)) zeros(n-1,1)]... 
  + [zeros(n-1,1) diag(1./delta(2:end))]...
  + [diag(1./delta(2:end-1),-1) zeros(n-1,1)];
end