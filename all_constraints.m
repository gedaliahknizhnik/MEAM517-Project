function [c, ceq, dC, dCeq] = all_constraints(z, N, nx, nu, dt, params)

% MODIFIED - Credit to Mat Halm
[ceq, dCeq] = dynamics_constraints(z, N, nx, nu, dt, params);

c = zeros(0,1);
dC = zeros(0,numel(z));

dC = sparse(dC);
dCeq = sparse(dCeq);

end