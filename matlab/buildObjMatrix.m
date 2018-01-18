%Build the objective matrix of the original problem
function J = buildObjMatrix(delta)
n = length(delta);
J = [diag(-1./delta(2:end)-1./delta(1:end-1)) zeros(n-1,1)]... 
  + [zeros(n-1,1) diag(1./delta(2:end))]...
  + [diag(1./delta(2:end-1),-1) zeros(n-1,1)];
end