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


%A small utility function to write latex code from a matrix
%\BIN and \BOUT, are shortcuts for \begin{bmatrix} and \end{bmatrix}
function s = matrix2latex(M)
[m,n] = size(M);
S = num2str(M(:,1));
delim = repmat(' & ',[m,1]);
for i=2:n
  S = [S delim num2str(M(:,i))];
end
l = size(S,2);
s = ['\BIN' repmat(' ',[1,l+1]);
     repmat(' ',[m,2]) S repmat(' \\',[m,1]); 
     '\BOUT' repmat(' ',[1,l])];
end