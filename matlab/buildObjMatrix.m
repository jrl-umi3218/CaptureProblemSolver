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