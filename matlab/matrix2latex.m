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